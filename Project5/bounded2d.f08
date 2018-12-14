program project5
  implicit none

    character(15) :: filename !HovmullerBoundedGaus.dat
    character(4) :: func_type
    character(3) :: hovmuller
    integer :: T,n,i,j,grid_size
    real :: delta_t,delta_x, t1, x, sigma,y
    DOUBLE PRECISION :: pi
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: psi(:),zeta(:),zeta_temp(:),zeta_use(:),zeta_store(:)

    T = 150
    grid_size = 40
    n = grid_size + 1
    func_type = 'sine'
    sigma = 0.2
    hovmuller = 'yes'

    ALLOCATE (psi(n*n),zeta(n*n),zeta_temp(n*n),zeta_use(n*n),zeta_store(n*n))

    delta_t = 0.01
    delta_x = 1.0/grid_size
    pi = 4.D0*DATAN(1.D0)

    ! Initializing psi and zeta for first timestep
    ! Also setting up grid vector x
    if (func_type == 'sine') then
      write(*,*) 'hello'
      DO i=1,n
        x = (i-1)*delta_x
        DO j=1,n
          y = (j-1)*delta_x
          ! Bounded domain conditions
          if (i == n) THEN
              psi((i-1)*n + j) = 0
              zeta((i-1)*n + j) = 0
          elseif (j == n) THEN
              psi((i-1)*n + j) = 0
              zeta((i-1)*n + j) = 0
          ! Initial function
          else
            psi((i-1)*n + j) = sin(pi*y)*sin(4.0*pi*x)
            zeta((i-1)*n + j) = -17.0*pi**2*sin(pi*y)*sin(4.0*pi*x)
          endif
        ENDDO
      ENDDO
    ! elseif (func_type == 'gaus') then
    !   write(*,*) 'nope'
    !   DO i=1,n
    !     x = (i-1)*delta_x
    !     psi(i) = exp(-((x-0.5)/sigma)**2)
    !     zeta(i) = (4.0*((x - 0.5)/(sigma*sigma))*((x - 0.5)/(sigma*sigma)) - 2.0/(sigma*sigma))&
    !     *exp(-((x-0.5)/sigma)*((x-0.5)/sigma))
    !   ENDDO
    endif


    ! Start time
    t1 = 0

    DO WHILE (t1 .LT. T)
      ! Use zeta to approximate first time/timestep
      if (t1 == 0) then
        DO i=2,n-1
          DO j = 2,n-1
            zeta_use((i-1)*n + j) = zeta((i-1)*n + j)
            zeta((i-1)*n + j) = zeta((i-1)*n + j) + delta_t/delta_x*(psi((i-2)*n + j) - psi(i*n + j))
            zeta_store((i-1)*n + j) = zeta((i-1)*n + j)
          ENDDO
        ENDDO
        zeta_temp = -zeta
        ! Update psi for first timestep
        CALL bounded2dsolver(zeta_temp,psi,n,delta_x)
      endif
      ! else
      !   ! Find zeta for next timestep
      !   CALL vorticity(zeta,zeta_use,zeta_store,psi,delta_x,delta_t,n)
      ! endif
      ! Insert zeta for next timestep to solve for psi at next timestep
      ! To get diff eq on form -u''(x) = f, have u''(x) = f in our case. Just putting in minus so dont have to change solver code.
      ! zeta_temp = -zeta*delta_x**2
      !CALL boundedsolver(zeta_temp,psi,n)
      ! For testing
      !IF (t1 == 0) THEN
      IF (t1 == 0) THEN
        11 format(A15)
        WRITE(filename,11)'wave2d_test.dat'

        open(unit=1,file=filename)
        write(1,*) psi
        close(1)
      ENDIF
      t1 = t1 + delta_t
    ENDDO

  DEALLOCATE (psi,zeta)
endprogram

SUBROUTINE bounded2dsolver(zeta, psi,n,delta_x)
  implicit none
  DOUBLE PRECISION, INTENT(INOUT) :: zeta(n*n),psi(n*n)
  INTEGER, INTENT(IN)             :: n
  INTEGER                         :: i,j,max_iteration,iteration
  REAL, INTENT(IN)                :: delta_x
  REAL                            :: delta_y,delta_x2,delta_y2,delta_x2y2,delta_x2_pluss_y2_2
  REAL                            :: diff,eps
  DOUBLE PRECISION                :: psi_store(n*n)

  delta_y = delta_x
  delta_x2 = delta_x**2
  delta_y2 = delta_y**2
  delta_x2y2 = delta_x**2*delta_x**2
  delta_x2_pluss_y2_2 = 2*(delta_x**2+delta_x**2)

  diff = 1E20
  eps = 1.0E-6
  max_iteration = 50
  iteration = 0

  DO WHILE (iteration .LE. max_iteration .AND. abs(diff) .GT. eps)
    diff = 0.0

    psi_store = psi
    DO i=2,n-1
      DO j = 2,n-1
        psi((i-1)*n + j) = (delta_y2*(psi_store((i)*n + j)+psi_store((i-2)*n + j)) + &
        delta_x2*(psi_store((i-1)*n + (j+1)) + psi_store((i)*n + (j-1))) - &
        delta_x2y2*zeta((i-1)*n + j))/delta_x2_pluss_y2_2
        diff = diff + (psi((i-1)*n + j) - psi_store((i-1)*n + j))
      ENDDO
    ENDDO
    iteration = iteration + 1
  ENDDO

end SUBROUTINE
