program project5
  implicit none

    character(25) :: filename !HovmullerBoundedGaus.dat
    character(4) :: func_type
    character(3) :: hovmuller
    integer :: T,n,i,j,grid_size
    real :: delta_t,delta_x, t1, x, sigma,y
    DOUBLE PRECISION :: pi
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: psi(:,:),zeta(:,:),zeta_temp(:,:),zeta_use(:,:),psi_store(:,:),zeta_store(:,:)

    T = 150
    grid_size = 40
    n = grid_size + 1
    func_type = 'sine'
    sigma = 0.2
    hovmuller = 'yes'

    ALLOCATE (psi(n,n),zeta(n,n),zeta_temp(n,n),zeta_use(n,n),psi_store(n,n),zeta_store(n,n))

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
          ! Initial function
          psi(i,j) = sin(4.0*pi*x)
          zeta(i,j) = -16.0*pi**2*sin(4.0*pi*x)
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
        CALL periodic2dsolver(zeta,psi,n,delta_x)
      else
        CALL vorticity(zeta,zeta_use,zeta_store,psi,delta_x,delta_t,n)
        CALL periodic2dsolver(zeta,psi,n,delta_x)
      endif

      ! For testing
      IF (t1 == 0) THEN
        11 format(A16)
        WRITE(filename,11)'2dperiodic_0.dat'

        open(unit=1,file=filename)
        write(1,*) psi
        close(1)
      ENDIF
      t1 = t1 + delta_t
    ENDDO

    12 format(A18)
    WRITE(filename,12)'2dperiodic_150.dat'

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

SUBROUTINE periodic2dsolver(zeta, psi,n,delta_x)
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

    DO WHILE ( iteration < max_iteration .AND. abs(diff) > eps)
      psi_store = psi

      !Assume dx == dy, so dx**2*dy**2 = dx**4, 2(dx**2 + dy**2) = 4dx**2

      !First we resolve the walls
      ! east west walls.
      DO j = 2, n-1
         psi(1,j) = (delta_y2*(psi_store(2,j)+psi_store(n,j)) + &
         delta_x2*(psi_store(1,j+1) + psi_store(1,j-1)) - &
         delta_x2y2*zeta(1,j))/delta_x2_pluss_y2_2

         psi(n,j) = (delta_y2*(psi_store(1,j)+psi_store(n-1,j)) + &
         delta_x2*(psi_store(n,j+1) + psi_store(n,j-1)) - &
         delta_x2y2*zeta(n,j))/delta_x2_pluss_y2_2
      END DO

      ! north south walls
      DO i = 2, n-1
        psi(i,1) = (delta_y2*(psi_store(i+1,1)+psi_store(i-1,1)) + &
        delta_x2*(psi_store(i,2) + psi_store(i,n)) - &
        delta_x2y2*zeta(i,1))/delta_x2_pluss_y2_2

        psi(i,n) = (delta_y2*(psi_store(i+1,n)+psi_store(i-1,n)) + &
        delta_x2*(psi_store(i,1) + psi_store(i,n-1)) - &
        delta_x2y2*zeta(i,n))/delta_x2_pluss_y2_2
      END DO

      !Second we resolve corners
      ! i = j = 1
      psi(1,1) = (delta_y2*(psi_store(2,1)+ psi_store(n,1)) + delta_x2*(psi_store(1,2) + psi_store(1,n)) - &
             delta_x2y2*zeta(1,1))/delta_x2_pluss_y2_2

      ! i = j = n
      psi(n,n) = (delta_y2*(psi_store(1,n)+ psi_store(n-1,n)) +  delta_x2*(psi_store(n,1) + psi_store(n,n-1)) - &
             delta_x2y2*zeta(n,n))/delta_x2_pluss_y2_2

      ! i = n,j = 1
      psi(n,1) = (delta_y2*(psi_store(1,1)+ psi_store(n-1,1)) + delta_x2*(psi_store(n,2) + psi_store(n,n)) - &
             delta_x2y2*zeta(n,1))/delta_x2_pluss_y2_2

      ! i = 1, j = n
      psi(1,n) = (delta_y2*(psi_store(2,n)+ psi_store(n,n)) + delta_x2*(psi_store(1,1) + psi_store(1,n-1)) - &
             delta_x2y2*zeta(1,n))/delta_x2_pluss_y2_2


      !Then the inner matrix
      DO i=2,n-1
        DO j = 2,n-1
          psi(i,j) = (delta_y2*(psi_store(i+1,j)+psi_store(i-1,j)) + &
          delta_x2*(psi_store(i,j+1) + psi_store(i,j-1)) - &
          delta_x2y2*zeta(i,j))/delta_x2_pluss_y2_2
        ENDDO
      ENDDO

      DO i=1,n
        DO j = 1,n
          diff = diff + (psi(i,j) - psi_store(i,j))
        ENDDO
      ENDDO
      iteration = iteration + 1
    END DO
  END SUBROUTINE
