program project5
  implicit none
    character(25) :: filename
    integer :: T,n,i,grid_size
    real :: delta_t,delta_x, t1, x
    DOUBLE PRECISION :: pi
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: psi(:),zeta(:),zeta_temp(:)

    T = 150
    grid_size = 40
    n = grid_size + 1

    ALLOCATE (psi(n),zeta(n),zeta_temp(n))

    delta_t = 0.01
    delta_x = 1.0/grid_size
    pi = 4.D0*DATAN(1.D0)

    ! Initializing psi and zeta for first timestep
    ! Also setting up grid vector x
    DO i=1,n
      x = (i-1)*delta_x
      psi(i) = sin(4.0*pi*x)
      zeta(i) = -(4*pi)**2*sin(4.0*pi*x)
    ENDDO
    ! Initial boundary conditions
    psi(1) = 0; psi(n) = 0; zeta(1) = 0; zeta(n) = 0

    ! Start time
    t1 = 0
    ! Implementation over time
    DO WHILE (t1 .LT. T)
      ! Find zeta for next timestep
      CALL vorticity(n,zeta,psi,delta_x,delta_t)
      ! Insert zeta for next timestep to solve for psi at next timestep
      ! To get diff eq on form -u''(x) = f, have u''(x) = f in our case. Just putting in minus so dont have to change solver code.
      zeta_temp = zeta*delta_x**2
      CALL periodicsolver(zeta_temp,psi,n)
      t1 = t1 + delta_t
    ENDDO

    ! Write endresult of wave at time T
    10 format(A14)
    WRITE(filename,10)'periodic1d.dat'

    open(unit=1,file=filename)
    write(1,*) psi
    close(1)

DEALLOCATE (psi,zeta)

end program

SUBROUTINE vorticity(n,zeta,psi,delta_x,delta_t)
  implicit NONE
  ! Input arguments
  integer,intent(in) :: n
  real,intent(in) :: delta_x,delta_t
  double precision,intent(inout) :: zeta(n),psi(n)
  ! Subroutine only arguments
  integer :: i
  real :: alpha

! Updated boundary conditions based on pariodic boundaries.
alpha = 0.5*delta_t/delta_x
zeta(1) = alpha*(psi(n) - psi(2)) + zeta(1)
zeta(n) = zeta(1)

! Updated zeta in time
DO i=2,n-1
  zeta(i) = zeta(i) + alpha*(psi(i-1) - psi(i+1))
ENDDO

END SUBROUTINE

SUBROUTINE periodicsolver(zeta_temp,psi,n)
use F90library
implicit NONE
! Inout arguments
double precision, intent(inout) :: zeta_temp(n), psi(n)
! Subroutine-only arguments
double precision                :: a(n,n)
integer                         :: n, i, j
integer                         :: indx(n)
double precision                :: d

do i = 1,n
  do j = 1,n
    if (i == j) then
      a(i,j) = -2.0
    elseif (i==(j+1)) then
      a(i,j) = 1.0
    elseif (i==(j-1)) then
      a(i,j) = 1.0
    else
      a(i,j) = 0
    endif
  end do
end do

a(1,n) = 1.0; a(n,1) = 1.0

call lu_decompose(a,n,indx,d)
call lu_linear_equation(a,n,indx,zeta_temp)
psi = zeta_temp
END SUBROUTINE
