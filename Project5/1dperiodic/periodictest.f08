program project5
  implicit none
    character(13) :: filename
    integer :: T,n,i,grid_size
    real :: delta_t,delta_x, t1, x
    DOUBLE PRECISION :: pi
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: psi(:),zeta(:),zeta_temp(:)

    T = 150
    grid_size = 40
    n = grid_size + 1

    ALLOCATE (psi(n),zeta(n),zeta_temp(n))

    delta_t = 0.01
    delta_x = 1.0/(grid_size+1)
    pi = 4.D0*DATAN(1.D0)

    ! Initializing psi and zeta for t = 0
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
      zeta_temp = -zeta*delta_x**2
      ! Calculate psi for next timestep
      CALL periodicsolver(zeta_temp,psi,n)
      t1 = t1 + delta_t
    ENDDO

    ! Write endresult of wave at time T
    10 format(A13)
    WRITE(filename,10)'wave_func.dat'

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

! Updated boundary conditions based on pariodic boundaries.
zeta(1) = delta_t/(2*delta_x)*(psi(n) - psi(2)) + zeta(1)
zeta(n) = zeta(1)! delta_t/(2*delta_x)*(psi(n-1) - psi(1)) + zeta(n)

! Updated zeta in time
DO i=2,n-1
  zeta(i) = zeta(i) + 0.5*delta_t/delta_x*(psi(i-1) - psi(i+1))
ENDDO

END SUBROUTINE

SUBROUTINE periodicsolver(f, v, n)
    USE f90library
    DOUBLE PRECISION, INTENT(INOUT) :: f(n), v(n)
    INTEGER, INTENT(IN)             :: n
    DOUBLE PRECISION                :: A(n,n), d
    INTEGER                         :: i, j, indx(n)

    DO i = 1,n
      DO j = 1,n
        IF (i == j) THEN
          !Fills the DIAGONAL
          A(i,j) = -2.d0
        ELSE IF (abs(i-j)==(1)) THEN
          !Fills BELOW and ABOVE (Distanced 1 from the DIAGONAL)
          A(i,j) = 1.d0
        ELSE
          A(i,j) = 0.d0
        END IF
      END DO
    END DO
    A(1,n) = 1.d0
    A(n,1) = 1.d0

    !Calls the LU-decomposition-algorithm, A is overwritten.
    CALL lu_decompose(A,n,indx,d)
    !Calls the solver, using f- as righthand-side. f is overwritten with v
    CALL lu_linear_equation(A,n,indx,f)
    v = f

  END SUBROUTINE
