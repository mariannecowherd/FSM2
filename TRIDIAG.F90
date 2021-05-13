!-----------------------------------------------------------------------
! Solve tridiagonal matrix equation
!-----------------------------------------------------------------------
subroutine TRIDIAG(Nvec,Nmax,a,b,c,r,x)

implicit none

integer, intent(in) :: &
  Nvec,              &! Vector length
  Nmax                ! Maximum vector length

real, intent(in) :: &
  a(Nmax),           &! Below-diagonal matrix elements
  b(Nmax),           &! Diagonal matrix elements
  c(Nmax),           &! Above-diagonal matrix elements
  r(Nmax)             ! Matrix equation rhs

real, intent(out) :: &
  x(Nmax)             ! Solution vector

integer :: n, nrhs, ldb,info
real, allocatable :: dl(:),d(:),du(:),bCp(:,:)

n = Nmax
nrhs = 1
ldb = Nmax

allocate(dl(n))
allocate(d(n))
allocate(du(n))
allocate(bCp(n,1))

dl(:) = a(:)
d(:) = b(:)
du(:) = c(:)
bCp(:,1) = r(:)

call sgtsv(Nvec,nrhs,dl,d,du,bCp,Nvec,info)
x(:) = bCp(:,1)
end subroutine TRIDIAG
