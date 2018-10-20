PROGRAM error_euler
  IMPLICIT None
  integer :: i,good_n
  DOUBLE PRECISION:: dist_error,last_dist_error

  last_dist_error = 0
  good_n = 0

  open(300,file="dist_error_euler.dat")
  do i=10,400
    CALL euler(i,dist_error)
    write(300,*) dist_error
    if (good_n > 0) then
      continue
    else
      ! For first loop run, store current_dist_error
      if (i == 10) then
        last_dist_error = dist_error
      ! For all other loop runs
      else
        ! If change in distance error larger than 0.1%, keep searching
        if (abs(dist_error-last_dist_error)/last_dist_error*100 > 0.1) then
          last_dist_error = dist_error
        ! If error below 0.1%, store loop run as good enough amount of integration points.
        else
          good_n = i
          write(*,*) good_n
        endif
      endif
    endif
  enddo
  close(300)


END PROGRAM error_euler

SUBROUTINE euler(n,dist_error)
  IMPLICIT NONE
  DOUBLE PRECISION :: FinalTime, h, pi, FourPi2,x_diff,y_diff
  DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:):: time, x, y, vx, vy
  DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:):: r
  INTEGER,intent(in) :: n
  DOUBLE PRECISION,intent(out)::dist_error
  INTEGER::i
  CHARACTER(19) ::filename
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: posx,posy,velx,vely,ax,ay,Fx,Fy

  FinalTime = 1.0
  h = FinalTime/REAL(n)
  pi = 4.D0*DATAN(1.D0)
  FourPi2 = 4*pi*pi
  ALLOCATE(time(n+1),x(n+1),y(n+1),vx(n+1),vy(n+1),r(n+1))
  time(1) = 0.0
  x(1) = 1.0
  y(1) = 0.0
  vx(1) = 0.0
  vy(1) = 2.0*pi
  r(1) = sqrt((x(1)**2)+(y(1)**2))

  DO i=1,n
      time(i+1) = time(i) + h
      x(i+1) = x(i) + (h*vx(i))
      y(i+1) = y(i) + (h*vy(i))
      vx(i+1) = vx(i) - h*FourPi2*x(i+1)/(r(i)*r(i)*r(i))
      vy(i+1) = vy(i) - h*FourPi2*y(i+1)/(r(i)*r(i)*r(i))
      r(i+1) = sqrt((x(i+1)**2)+(y(i+1)**2))
  ENDDO
  x_diff = x(1)-x(n)
  y_diff = y(1)-y(n)
  dist_error = sqrt((x_diff**2)+(y_diff**2))
  DEALLOCATE(time,x,y,vx,vy,r)
END SUBROUTINE euler
