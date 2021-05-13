!-----------------------------------------------------------------------
! Solve matrix equation Ax = b for x by LU decomposition
!-----------------------------------------------------------------------
subroutine LUDCMP(N,A,b,x)
implicit none
integer, intent(in) :: &
  N                   ! Number of equations to solve

real, intent(in) :: &
 A(N,N),            & ! Matrix
 b(N)                 ! RHS of matrix equation

real, intent(out) :: &
 x(N)                 ! Solution of matrix equation

integer ::  nrhs, lda, ldb, info
real, allocatable  :: aCp(:,:), bCp(:,:)
integer,allocatable :: ipiv(:)
nrhs = 1
lda = N
ldb = N

allocate(aCp(N,N))
allocate(bCp(N,nrhs))
allocate(ipiv(N))

aCP(:,:) = A(:,:)
bCp(:,nrhs) = b(:)

call dgesv(N,nrhs,aCp,lda,ipiv,bCp,ldb,info) 
x(:) = bCp(:,1)
end subroutine LUDCMP
