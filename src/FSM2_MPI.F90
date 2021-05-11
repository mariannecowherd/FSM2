!----------------------------------------------------------------------!
! Flexible Snow Model (FSM version 2.1.0)                              !
!                                                                      !
! Richard Essery                                                       !
! School of GeoSciences                                                !
! University of Edinburgh                                              !
!----------------------------------------------------------------------!
program FSM2

#include "OPTS.h"

use mpi

#if PROFNC == 1
use netcdf
#endif

use CONSTANTS, only: &
  hcon_air,          &! Thermal conductivity of air (W/m/K)
  hcon_clay,         &! Thermal conductivity of clay (W/m/K)
  hcon_sand           ! Thermal conductivity of sand (W/m/K)

use IOUNITS, only: &
  ucan,              &! Subcanopy diagnostics file unit number
  udmp,              &! Start / dump file unit number
  uflx,              &! Flux output file unit number
  umet,              &! Meteorological driving file unit number
  usta                ! State output file unit number

use LAYERS, only: &
  Dzsnow,            &! Minimum snow layer thicknesses (m)
  Dzsoil,            &! Soil layer thicknesses (m)
  fvg1,              &! Fraction of vegetation in upper canopy layer
  Ncnpy,             &! Number of canopy layers
  Nsmax,             &! Maximum number of snow layers
  Nsoil,             &! Number of soil layers
  zsub                ! Subcanopy wind speed diagnostic height (m)

use PARAMETERS, only: &
  fcly,              &! Soil clay fraction
  fsnd,              &! Soil sand fraction
  rgr0                ! Fresh snow grain radius (m)

use SOILPROPS, only: &
  b,                 &! Clapp-Hornberger exponent
  hcap_soil,         &! Volumetric heat capacity of dry soil (J/K/m^3)
  hcon_soil,         &! Thermal conductivity of dry soil (W/m/K)
  sathh,             &! Saturated soil water pressure (m)
  Vcrit,             &! Volumetric soil moisture at critical point
  Vsat                ! Volumetric soil moisture at saturation

implicit none

! Grid dimensions
integer :: &
  Ncols,             &! Number of columns in grid
  Nrows               ! Number of rows in grid

! Site characteristics
real :: &
  lat,               &! Latitude (radians)
  noon                ! Time of solar noon (hour)

! Meteorological driving data
character(len=70) :: &
  met_file            ! Meteorological driving file name
integer :: &
  year,              &! Year
  month,             &! Month of year
  day                 ! Day of month
logical :: EoF        ! End-of-file flag
real :: &
  dt,                &! Timestep (s)
  elev,              &! Solar elevation (radians)
  hour,              &! Hour of day
  zT,                &! Temperature and humidity measurement height (m)
  zU                  ! Wind speed measurement height (m)
real, allocatable :: &
  LW(:,:),           &! Incoming longwave radiation (W/m^2)
  Ps(:,:),           &! Surface pressure (Pa)
  Qa(:,:),           &! Specific humidity (kg/kg)
  Rf(:,:),           &! Rainfall rate (kg/m^2/s)
  Sdif(:,:),         &! Diffuse shortwave radiation (W/m^2)
  Sdir(:,:),         &! Direct-beam shortwave radiation (W/m^2)
  Sf(:,:),           &! Snowfall rate (kg/m^2/s)
  Ta(:,:),           &! Air temperature (K)
  trans(:,:),        &! Wind-blown snow transport rate (kg/m^2/s)
  Ua(:,:)             ! Wind speed (m/s)

! Model state variables
integer, allocatable :: &
  Nsnow(:,:)          ! Number of snow layers
real, allocatable :: &
  albs(:,:),         &! Snow albedo
  Tsrf(:,:),         &! Snow/ground surface temperature (K)
  Dsnw(:,:,:),       &! Snow layer thicknesses (m)
  Qcan(:,:,:),       &! Canopy air space humidities
  Rgrn(:,:,:),       &! Snow layer grain radii (m)
  Sice(:,:,:),       &! Ice content of snow layers (kg/m^2)
  Sliq(:,:,:),       &! Liquid content of snow layers (kg/m^2)
  Sveg(:,:,:),       &! Snow mass on vegetation layers (kg/m^2)
  Tcan(:,:,:),       &! Canopy air space temperatures (K)
  Tsnow(:,:,:),      &! Snow layer temperatures (K)
  Tsoil(:,:,:),      &! Soil layer temperatures (K)
  Tveg(:,:,:),       &! Vegetation layer temperatures (K)
  Vsmc(:,:,:),       &! Volumetric moisture content of soil layers
  fsat(:),           &! Initial soil layer moisture/saturation
  Tprf(:)             ! Initial soil layer temperatures (K)

! Diagnostics
real, allocatable :: &
  H(:,:),            &! Sensible heat flux to the atmosphere (W/m^2)
  LE(:,:),           &! Latent heat flux to the atmosphere (W/m^2)
  LWout(:,:),        &! Outgoing LW radiation (W/m^2)
  LWsub(:,:),        &! Subcanopy downward LW radiation (W/m^2)
  Melt(:,:),         &! Surface melt rate (kg/m^2/s)
  Roff(:,:),         &! Runoff from snow (kg/m^2/s)
  snd(:,:),          &! Snow depth (m)
  snw(:,:),          &! Total snow mass on ground (kg/m^2)
  subl(:,:),         &! Sublimation rate (kg/m^2/s)
  svg(:,:),          &! Total snow mass on vegetation (kg/m^2)
  SWout(:,:),        &! Outgoing SW radiation (W/m^2)
  SWsub(:,:),        &! Subcanopy downward SW radiation (W/m^2)
  Usub(:,:),         &! Subcanopy wind speed (m/s)
  Wflx(:,:,:)         ! Water flux into snow layer (kg/m^2/s)

! Vegetation characteristics
character(len=70) :: &
  alb0_file,         &! Snow-free ground albedo map file name
  fsky_file,         &! Skyview fraction map file name
  vegh_file,         &! Canopy height map file name
  VAI_file            ! Vegetation area index map file name
real, allocatable :: &
  alb0(:,:),         &! Snow-free ground albedo
  fsky(:,:),         &! Skyview not obstructed by remote vegetation
  vegh(:,:),         &! Canopy height (m)
  VAI(:,:)            ! Vegetation area index

! Start and dump file names
character(len=70) :: &
  dump_file,         &! End dump file name
  runid,             &! Run identifier
  start_file          ! Start file name

! NetCDF variables
integer :: &
  ncid,              &! Dataset ID
  rec,               &! Record number
  status,            &! Error status
  varid(17)           ! Variable IDs

! Counters
integer :: &
  i,                 &! Grid row counter
  j,                 &! Grid column counter
  k                   ! Soil layer counter


  integer ( kind = 4 ) error
  integer ( kind = 4 ) id
  integer ( kind = 4 ) p
  real ( kind = 8 ) wtime
  character(100):: id_string
  integer :: rows_per
  integer :: my_start
  integer :: my_end




namelist    /drive/ met_file,dt,lat,noon,zT,zU
namelist /gridpnts/ Ncols,Nrows,Nsmax,Nsoil
namelist /gridlevs/ Dzsnow,Dzsoil,fvg1,zsub
namelist  /initial/ fsat,Tprf,start_file
namelist  /outputs/ dump_file,runid
namelist      /veg/ alb0,fsky,vegh,VAI,  &
                    alb0_file,fsky_file,vegh_file,VAI_file

#if SETPAR == 1
call FSM2_PARAMS
#endif



!
!  Initialize MPI.
!
  call MPI_Init ( error )
!
!  Get the number of processes.
!
  call MPI_Comm_size ( MPI_COMM_WORLD, p, error )
!
!  Get the individual process ID.
!
  call MPI_Comm_rank ( MPI_COMM_WORLD, id, error )

    wtime = MPI_Wtime ( )

! Grid dimensions
Ncols = 20
Nrows = 20
Nsmax = 3
Nsoil = 4
! read(5,gridpnts)

rows_per = (Nrows / p)
my_start = rows_per * id + 1
my_end = rows_per * (id+1) + 1
if (my_start > Nrows ) stop
if (my_end > Nrows ) then
  my_end = Nrows
end if


! Canopy, snow and soil layers
#if CANMOD == 1
Ncnpy = 1
#endif
#if CANMOD == 2
Ncnpy = 2
#endif
fvg1 = 0.5
zsub = 1.5
allocate(Dzsnow(Nsmax))
allocate(Dzsoil(Nsoil))
if (Nsmax == 3) Dzsnow = (/0.1, 0.2, 0.4/)
if (Nsoil == 4) Dzsoil = (/0.1, 0.2, 0.4, 0.8/)
! read(5, gridlevs)

! Site and driving data characteristics
met_file = 'met'
dt = 3600
lat = 0
noon = 12
zT = 2
zU = 10
! read(5,drive)

met_file = 'met_Sod_1314.txt'
zT = 18
zU = 18
lat = 67.37
noon = 10

open(umet, file = met_file)
read(umet,*) year,month,day,hour
rewind(umet)
lat = (3.14159/180)*lat

! Allocate driving data arrays
allocate(LW(Ncols,rows_per))
allocate(Ps(Ncols,rows_per))
allocate(Qa(Ncols,rows_per))
allocate(Rf(Ncols,rows_per))
allocate(Sdif(Ncols,rows_per))
allocate(Sdir(Ncols,rows_per))
allocate(Sf(Ncols,rows_per))
allocate(Ta(Ncols,rows_per))
allocate(trans(Ncols,rows_per))
allocate(Ua(Ncols,rows_per))
trans(:,:) = 0

! Vegetation characteristics from defaults, namelist or named map files
allocate(alb0(Ncols,rows_per))
allocate(fsky(Ncols,rows_per))
allocate(vegh(Ncols,rows_per))
allocate(VAI(Ncols,rows_per))
alb0_file = 'none'
fsky_file = 'none'
vegh_file = 'none'
VAI_file  = 'none'
alb0 = 0.2
fsky = 1
vegh = 0
VAI  = 0
! read(5,veg)
vegh = 15
VAI =  1
if (alb0_file /= 'none') call FSM2_MAP(alb0_file,Ncols,rows_per,alb0)
if (fsky_file /= 'none') call FSM2_MAP(fsky_file,Ncols,rows_per,fsky)
if (vegh_file /= 'none') call FSM2_MAP(vegh_file,Ncols,rows_per,vegh)
if (VAI_file  /= 'none') call FSM2_MAP(VAI_file,Ncols,rows_per,VAI)

! Soil properties
b = 3.1 + 15.7*fcly - 0.3*fsnd
hcap_soil = (2.128*fcly + 2.385*fsnd)*1e6 / (fcly + fsnd)
sathh = 10**(0.17 - 0.63*fcly - 1.58*fsnd)
Vsat = 0.505 - 0.037*fcly - 0.142*fsnd
Vcrit = Vsat*(sathh/3.364)**(1/b)
hcon_soil = (hcon_air**Vsat) * ((hcon_clay**fcly)*(hcon_sand**(1 - fcly))**(1 - Vsat))

! Allocate state variable arrays
allocate(albs(Ncols,rows_per))
allocate(Nsnow(Ncols,rows_per))
allocate(Tsrf(Ncols,rows_per))
allocate(Dsnw(Nsmax,Ncols,rows_per))
allocate(Qcan(Ncnpy,Ncols,rows_per))
allocate(Rgrn(Nsmax,Ncols,rows_per))
allocate(Sice(Nsmax,Ncols,rows_per))
allocate(Sliq(Nsmax,Ncols,rows_per))
allocate(Sveg(Ncnpy,Ncols,rows_per))
allocate(Tcan(Ncnpy,Ncols,rows_per))
allocate(Tsnow(Nsmax,Ncols,rows_per))
allocate(Tsoil(Nsoil,Ncols,rows_per))
allocate(Tveg(Ncnpy,Ncols,rows_per))
allocate(Vsmc(Nsoil,Ncols,rows_per))

! Default initialization of state variables
albs  = 0.8
Dsnw  = 0
Nsnow = 0
Qcan  = 0
Rgrn  = rgr0
Sice  = 0
Sliq  = 0
Sveg  = 0
Tcan  = 285
Tsnow = 273
Tveg  = 285
! Missing values for vegetation at non-forest points
do k = 1, Ncnpy
  where(VAI==0) Sveg(k,:,:) = -999./Ncnpy
  where(VAI==0) Tveg(k,:,:) = -999
end do

! Initial soil profiles from namelist
allocate(fsat(Nsoil))
allocate(Tprf(Nsoil))
fsat = 0.5
Tprf = 285
start_file = 'Sod_1314_dump'
! read(5,initial)
do k = 1, Nsoil
  Tsoil(k,:,:) = Tprf(k)
  Vsmc(k,:,:) = fsat(k)*Vsat
end do
Tsrf = Tsoil(1,:,:)


! Initialize state variables from a named start file
if (start_file /= 'none') then
  open(udmp,file = start_file)
  read(udmp,*) albs
  read(udmp,*) Dsnw
  read(udmp,*) Nsnow
  read(udmp,*) Qcan
  read(udmp,*) Rgrn
  read(udmp,*) Sice
  read(udmp,*) Sliq
  read(udmp,*) Sveg
  read(udmp,*) Tcan
  read(udmp,*) Tsnow
  read(udmp,*) Tsoil
  read(udmp,*) Tsrf
  read(udmp,*) Tveg
  read(udmp,*) Vsmc
  close(udmp)
end if

! Allocate diagnostic output arrays
allocate(H(Ncols,rows_per))
allocate(LE(Ncols,rows_per))
allocate(LWout(Ncols,rows_per))
allocate(LWsub(Ncols,rows_per))
allocate(Melt(Ncols,rows_per))
allocate(Roff(Ncols,rows_per))
allocate(snd(Ncols,rows_per))
allocate(snw(Ncols,rows_per))
allocate(subl(Ncols,rows_per))
allocate(svg(Ncols,rows_per))
allocate(SWout(Ncols,rows_per))
allocate(SWsub(Ncols,rows_per))
allocate(Usub(Ncols,rows_per))
allocate(Wflx(Nsmax,Ncols,rows_per))

! Output files
dump_file = 'dump'
runid = 'none'
! read(5,outputs)
runid = 'five_'
if (runid == 'none') runid = ''
#if PROFNC == 1
if (Ncols*Nrows>1) stop 'NetCDF output only available for Nrows = Ncols = 1'
call FSM2_PREPNC(runid,year,month,day,hour,ncid,rec,varid)
#else
write (id_string,'(I0)') id
if (maxval(VAI) > 0) open(ucan, file = 'out/subc'//trim(id_string)//'.txt')
open(uflx, file = 'out/flux'//trim(id_string)//'.txt')
open(usta, file = 'out/stat'//trim(id_string)//'.txt')
#endif



! Run the model

EoF = .false.


do
  call FSM2_DRIVE(Ncols,rows_per,fsky,lat,noon,                           &
                  year,month,day,hour,elev,EoF,                        &
                  LW,Ps,Qa,Rf,Sdif,Sdir,Sf,Ta,Ua)
  if (EoF) goto 1

  do i = 1, rows_per
  do j = 1, Ncols
    call FSM2_TIMESTEP(                                                &
                       ! Driving variables                             &
                       dt,elev,zT,zU,LW(j,i),Ps(j,i),Qa(j,i),          &
                       Rf(j,i),Sdif(j,i),Sdir(j,i),Sf(j,i),            &
                       Ta(j,i),trans(j,i),Ua(j,i),                     &
                       ! Vegetation characteristics                    &
                       alb0(j,i),vegh(j,i),VAI(j,i),                   &
                       ! State variables                               &
                       albs(j,i),Tsrf(j,i),Dsnw(:,j,i),Nsnow(j,i),     &
                       Qcan(:,j,i),Rgrn(:,j,i),Sice(:,j,i),            &
                       Sliq(:,j,i),Sveg(:,j,i),Tcan(:,j,i),            &
                       Tsnow(:,j,i),Tsoil(:,j,i),Tveg(:,j,i),          &
                       Vsmc(:,j,i),                                    &
                       ! Diagnostics                                   &
                       H(j,i),LE(j,i),LWout(j,i),LWsub(j,i),           &
                       Melt(j,i),Roff(j,i),snd(j,i),snw(j,i),          &
                       subl(j,i),svg(j,i),SWout(j,i),SWsub(j,i),       &
                       Usub(j,i),Wflx(:,j,i)                           )
  end do
  end do


  !  call FSM2_OUTPUT(Ncols,rows_per,year,month,day,hour,                    &
  !                   H,LE,LWout,LWsub,Melt,Roff,snd,snw,subl,svg,SWout,  &
  !                   SWsub,Tsoil,Tsrf,Tveg,Usub,VAI)


end do
1 continue
! Write out state variables at end of run
open(udmp,file = 'out/'//trim(dump_file)//trim(id_string))

write(udmp,*) albs
write(udmp,*) Dsnw
write(udmp,*) Nsnow
write(udmp,*) Qcan
write(udmp,*) Rgrn
write(udmp,*) Sice
write(udmp,*) Sliq
write(udmp,*) Sveg
write(udmp,*) Tcan
write(udmp,*) Tsnow
write(udmp,*) Tsoil
write(udmp,*) Tsrf
write(udmp,*) Tveg
write(udmp,*) Vsmc
close(udmp)





    wtime = MPI_Wtime ( ) - wtime
    write ( *, '(a)' ) ''
    write ( *, '(a,i1,2x,a,g14.6,a)' ) &
      'P', id, '  Elapsed wall clock time = ', wtime, ' seconds.'
    write (*, * ) 'P', id, ' started:', my_start
    write (*, *) 'ended', my_end

    call MPI_Finalize ( error )

end program FSM2
