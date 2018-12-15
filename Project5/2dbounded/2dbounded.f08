program project5
  implicit none

    character(25) :: filename !HovmullerBoundedGaus.dat
    character(4) :: func_type
    character(3) :: hovmuller
    integer :: T,n,i,j,grid_size
    real :: delta_t,delta_x, t1, x, sigma,y
    DOUBLE PRECISION :: pi
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: psi(:,:),zeta(:,:),zeta_temp(:,:),zeta_use(:,:),zeta_store(:,:)

    T = 150
    grid_size = 40
    n = grid_size + 1
    func_type = 'sine'
    hovmuller = 'yes'

    ALLOCATE (psi(n,n),zeta(n,n),zeta_temp(n,n),zeta_use(n,n),zeta_store(n,n))

    delta_t = 0.025
    delta_x = 1.0/grid_size
    pi = 4.D0*DATAN(1.D0)

    ! Initializing psi and zeta for first timestep
    ! Also setting up grid vector x
    if (func_type == 'sine') then
      DO i=1,n
        x = (i-1)*delta_x
        DO j=1,n
          y = (j-1)*delta_x
          ! Bounded domain conditions
          if (i == n) THEN
              psi(i,j) = 0
              zeta(i,j) = 0
          elseif (j == n) THEN
              psi(i,j) = 0
              zeta(i,j) = 0
          elseif (i == 0) THEN
              psi(i,j) = 0
              zeta(i,j) = 0
          elseif (j == 0) THEN
              psi(i,j) = 0
              zeta(i,j) = 0
          ! Initial function
          else
            psi(i,j) = sin(pi*y)*sin(4.0*pi*x)
            zeta(i,j) = -(4.0*pi)**2*sin(pi*y)*sin(4.0*pi*x)
          endif
        ENDDO
      ENDDO
    endif

    ! Start time
    t1 = 0

    DO WHILE (t1 .LT. T)
      ! Use zeta to approximate first time/timestep
      if (t1 == 0) then
        zeta_use = zeta
        DO i=2,n-1
          DO j = 2,n-1
            zeta(i,j) = zeta_use(i,j) + 0.5*delta_t/delta_x*(psi(i-1,j) - psi(i+1,j))
          ENDDO
        ENDDO
        zeta_store = zeta
        ! Update psi for first timestep
        CALL bounded2dsolver(zeta,psi,n,delta_x)
      else
        CALL vorticity(zeta,zeta_use,zeta_store,psi,delta_x,delta_t,n)

        CALL bounded2dsolver(zeta,psi,n,delta_x)
      endif

      IF (t1 == 0) THEN
        11 format(A15)
        WRITE(filename,11)'2dbounded_0.dat'

        open(unit=1,file=filename)
        write(1,*) psi
        close(1)
      ENDIF
      IF (abs(t1 - 50) < 0.001) THEN
        12 format(A16)
        WRITE(filename,12)'2dbounded_50.dat'

        open(unit=1,file=filename)
        write(1,*) psi
        close(1)
      ENDIF

      IF (abs(t1 - 100) < 0.001) THEN
        13 format(A17)
        WRITE(filename,13)'2dbounded_100.dat'

        open(unit=1,file=filename)
        write(1,*) psi
        close(1)
      ENDIF
      t1 = t1 + delta_t
    ENDDO

    14 format(A17)
    WRITE(filename,14)'2dbounded_150.dat'

    open(unit=1,file=filename)
    write(1,*) psi
    close(1)

  DEALLOCATE (psi,zeta)
endprogram

SUBROUTINE vorticity(zeta,zeta_use,zeta_store,psi,delta_x,delta_t,n)
  implicit NONE
  ! Input arguments
  integer,intent(in) :: n
  real,intent(in) :: delta_x,delta_t
  double precision, dimension(n,n), intent(inout) :: zeta,psi,zeta_use,zeta_store
  ! Subroutine only arguments
  integer :: i,j

  DO i=2,n-1
    DO j = 2,n-1
      zeta(i,j) = zeta_use(i,j) + delta_t/delta_x*(psi(i-1,j) - psi(i+1,j))
    ENDDO
  ENDDO
  ! Update variables
  zeta_use = zeta_store
  zeta_store = zeta

END SUBROUTINE

SUBROUTINE bounded2dsolver(zeta, psi,n,delta_x)
  implicit none
  DOUBLE PRECISION, INTENT(INOUT) :: zeta(n,n),psi(n,n)
  INTEGER, INTENT(IN)             :: n
  INTEGER                         :: i,j,max_iteration,iteration
  REAL, INTENT(IN)                :: delta_x
  REAL                            :: delta_y,delta_x2,delta_y2,delta_x2y2,delta_x2_pluss_y2_2
  REAL                            :: diff,eps
  DOUBLE PRECISION                :: psi_store(n,n)

  delta_y = delta_x
  delta_x2 = delta_x**2
  delta_y2 = delta_x**2
  delta_x2y2 = delta_x**4
  delta_x2_pluss_y2_2 = 4.0*delta_x2

  diff = 100.0
  eps = 1.0E-6
  max_iteration = 100
  iteration = 0

  DO WHILE (iteration .LE. max_iteration .AND. abs(diff) .GT. eps)
    diff = 0.0

    psi_store = psi
    DO i=2,n-1
      DO j = 2,n-1
        psi(i,j) = (delta_y2*(psi_store(i+1,j)+psi_store(i-1,j)) + &
        delta_x2*(psi_store(i,j+1) + psi_store(i,j-1)) - &
        delta_x2y2*zeta(i,j))/delta_x2_pluss_y2_2
        diff = diff + (psi(i,j) - psi_store(i,j))
      ENDDO
    ENDDO
    iteration = iteration + 1
  ENDDO

end SUBROUTINE
