program project5
    ! Program for solving the wave equation in the bounded domain by using the leapfrog time scheme.
    ! Subroutine vorticity advances the vorticity to the next timestep
    ! Subroutine boundedsolver uses the advanced vorticity to calculate the next psi.

    ! Can also chose between plotting a sine and gausian wave by changing the variable func_type
    ! And if variable hovmuller = 'yes' one gets a datafile for plotting hovmuller diagrams.

    implicit none
    character(25) :: filename !HovmullerBoundedGaus.dat
    character(4) :: func_type
    character(3) :: hovmuller
    integer :: T,n,i,grid_size
    real :: delta_t,delta_x, t1, x, sigma
    DOUBLE PRECISION :: pi
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: psi(:),zeta(:),zeta_temp(:),zeta_use(:),zeta_store(:)

    T = 150
    grid_size = 40
    n = grid_size + 1
    func_type = 'gaus'
    sigma = 0.2
    !sigma = 0.1
    hovmuller = 'yes'

    ALLOCATE (psi(n),zeta(n),zeta_temp(n),zeta_use(n),zeta_store(n))

    delta_t = 0.02
    delta_x = 1.0/grid_size
    pi = 4.D0*DATAN(1.D0)

    ! Initializing psi and zeta for first timestep
    ! Also setting up grid vector x
    if (func_type == 'sine') then
      DO i=1,n
        x = (i-1)*delta_x
        psi(i) = sin(4.d0*pi*x)
        zeta(i) = -(4.d0*pi)**2*sin(4.d0*pi*x)
      ENDDO
    elseif (func_type == 'gaus') then
      DO i=1,n
        x = (i-1)*delta_x
        psi(i) = exp(-((x-0.5)/sigma)**2)
        zeta(i) = (4.0*((x - 0.5)/(sigma*sigma))*((x - 0.5)/(sigma*sigma)) - 2.0/(sigma*sigma))&
        *exp(-((x-0.5)/sigma)*((x-0.5)/sigma))
      ENDDO
    endif
    ! Initial boundary conditions
    psi(1) = 0; psi(n) = 0; zeta(1) = 0; zeta(n) = 0

    ! Start time
    t1 = 0
    ! Implementation over time
    ! Write endresult of wave at time T
    if (func_type == 'sine') then
      13 format(A24)
      WRITE(filename,13)'HovmullerBoundedSine.dat'
    else if (func_type == 'gaus' .AND. sigma == 0.2) then
      14 format(A25)
      WRITE(filename,14)'HovmullerBoundedGaus2.dat'
    else if (func_type == 'gaus' .AND. sigma == 0.1) then
      15 format(A25)
      WRITE(filename,15)'HovmullerBoundedGaus1.dat'
    endif
    open(unit=5,file=filename)
    DO WHILE (t1 .LT. T)
      if (hovmuller == 'yes') then
        write(5,*) psi
      endif
      ! Use Euler forward as first approximation
      if (t1 == 0) then
        DO i=2,n-1
          ! Euler forward to get enough zeta's for leapfrog
          zeta_use(i) = zeta(i)
          zeta(i) = zeta(i) + delta_t/delta_x*(psi(i-1) - psi(i+1))
          zeta_store(i) = zeta(i)
        ENDDO
        zeta_temp = -zeta*delta_x**2
        CALL boundedsolver(zeta_temp,psi,n)
      else
        ! Find zeta for next timestep
        CALL vorticity(zeta,zeta_use,zeta_store,psi,delta_x,delta_t,n)
        ! Insert zeta for next timestep to solve for psi at next timestep
        ! To get diff eq on form -u''(x) = f, have u''(x) = f in our case. Just putting in minus so dont have to change solver code.
        zeta_temp = -zeta*delta_x**2
        CALL boundedsolver(zeta_temp,psi,n)
      endif
      t1 = t1 + delta_t
    ENDDO
    close(5)

    ! Write endresult of wave at time T
    10 format(A18)
    WRITE(filename,10)'1dbounded_leap.dat'

    open(unit=1,file=filename)
    write(1,*) psi
    close(1)

DEALLOCATE (psi,zeta)
end program

SUBROUTINE vorticity(zeta,zeta_use,zeta_store,psi,delta_x,delta_t,n)
  implicit NONE
  ! Input arguments
  integer,intent(in) :: n
  real,intent(in) :: delta_x,delta_t
  double precision, dimension(n), intent(inout) :: zeta,psi,zeta_use,zeta_store
  ! Subroutine only arguments
  integer :: i

! Updated boundary conditions based on bounded boundaries.
zeta(1) = 0; zeta(n) = 0

! Updated zeta
DO i=2,n-1
  zeta(i) = zeta_use(i) + delta_t/delta_x*(psi(i-1) - psi(i+1))
ENDDO

! Update variables
zeta_use = zeta_store
zeta_store = zeta

END SUBROUTINE

SUBROUTINE boundedsolver(zeta_temp,psi,n) ! tridiag_solv()
  implicit none
  DOUBLE PRECISION, INTENT(INOUT) :: zeta_temp(n),psi(n)
  DOUBLE PRECISION                :: b(n)!,a(n),frac(n)
  INTEGER, INTENT(IN)             :: n
  INTEGER                         :: i

  b(1)=2; b(n)=2; psi(1)=0; psi(n)=0
! Algorithm
    ! ! Setting up b
    DO i=2,n-1
      b(i) = (i)/(i-1.0) ! 2 flops
    ENDDO
    ! Forward substitution
    DO i=3,n-1
      zeta_temp(i) = zeta_temp(i) + zeta_temp(i-1)/b(i-1)! 2 flops
    ENDDO

    ! Backward substitution
    psi(n-1) = zeta_temp(n-1)/b(n-1)

    DO i=n-1,2,-1
      psi(i) = (zeta_temp(i)+psi(i+1))/b(i) ! 2 flops
    ENDDO
  ! End Algorithm
END SUBROUTINE
