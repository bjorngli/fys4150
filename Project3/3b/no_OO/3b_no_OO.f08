PROGRAM calculate

  CALL euler()
  CALL verlet()

END PROGRAM calculate

SUBROUTINE euler()
  IMPLICIT NONE
  DOUBLE PRECISION :: FinalTime, h, pi, FourPi2, hFourPi2
  DOUBLE PRECISION, ALLOCATABLE :: x(:), y(:), vx(:), vy(:), r(:)
  INTEGER :: n, i
  REAL :: start, finish
  CHARACTER(19) ::filename

  FinalTime = 1000.
  n = 100000
  h = FinalTime/n
  pi = 4.D0*DATAN(1.D0)
  FourPi2 = 4*pi*pi
  hFourPi2 = h*FourPi2
  ALLOCATE(x(n+1),y(n+1),vx(n+1),vy(n+1),r(n+1))
  x(1) = 1.0
  y(1) = 0.0
  vx(1) = 0.0
  vy(1) = 2.0*pi
  r(1) = sqrt((x(1)**2)+(y(1)**2))
  call cpu_time(start)
  DO i=1,n+1
      x(i+1) = x(i) + (h*vx(i))
      y(i+1) = y(i) + (h*vy(i))
      vx(i+1) = vx(i) - hFourPi2*x(i+1)/(r(i)**2)
      vy(i+1) = vy(i) - hFourPi2*y(i+1)/(r(i)**2)
      r(i+1) = sqrt((x(i+1)**2)+(y(i+1)**2))
  ENDDO
  call cpu_time(finish)
  WRITE(*,*) 'euler: ', finish-start

  10 format(A19)
  WRITE(filename,10)'xpos_euler_noOO.dat'
  OPEN(unit=1,file=filename)
  WRITE(1,*) x
  CLOSE(1)

  11 format(A19)
  WRITE(filename,11)'ypos_euler_noOO.dat'
  OPEN(unit=1,file=filename)
  WRITE(1,*) y
  CLOSE(1)

  DEALLOCATE(x,y,vx,vy,r)
END SUBROUTINE euler

SUBROUTINE verlet()
  IMPLICIT NONE
  DOUBLE PRECISION,ALLOCATABLE :: x(:),y(:),vx(:),vy(:),r(:)
  INTEGER :: n,i ! Number of planets
  REAL :: start1, finish1
  DOUBLE PRECISION :: solar_mass, earth_mass,pi,FourPi2,FinalTime,h,solar_mass_value
  DOUBLE PRECISION :: ri2, hsqrd2, hdiv2, FourPi2x, FourPi2y
  CHARACTER(20) ::filename

  ! Variables
  FinalTime = 1000.0
  n = 100000
  h = FinalTime/n
  pi = 4.D0*DATAN(1.D0)
  FourPi2 = 4*pi*pi
  solar_mass_value = 1.989e30
  solar_mass = 1 ! 1 AU
  earth_mass = 6e24/solar_mass_value ! earth mass in units of AU

  ALLOCATE(x(n+2),y(n+2),vx(n+2),vy(n+2),r(n+2))

  ! Initial conditions
  hsqrd2=h**2/2
  hdiv2=h/2
  x(1) = 1.0
  y(1) = 0.0
  vx(1) = 0.0
  vy(1) = 2.0*pi
  r(1) = sqrt((x(1)**2)+(y(1)**2))
  call cpu_time(start1)
  DO i=1,n+1
      ri2=r(i)**2
      FourPi2x=FourPi2*x(i)
      FourPi2y=FourPi2*y(i)
      x(i+1) = x(i) + h*vx(i) - hsqrd2*FourPi2x/(ri2)
      y(i+1) = y(i) + h*vy(i) - hsqrd2*FourPi2y/(ri2)
      r(i+1) = sqrt((x(i+1)**2)+(y(i+1)**2))
      vx(i+1) = vx(i) - hdiv2*((FourPi2*x(i+1))/(r(i+1)**2)+FourPi2x/(ri2))
      vy(i+1) = vy(i) - hdiv2*((FourPi2*y(i+1))/(r(i+1)**2)+FourPi2y/(ri2))
  ENDDO
  call cpu_time(finish1)
  WRITE(*,*) 'verlet: ', finish1-start1

  10 format(A20)
  WRITE(filename,10)'xpos_verlet_noOO.dat'
  OPEN(unit=1,file=filename)
  WRITE(1,*) x
  CLOSE(1)

  11 format(A20)
  WRITE(filename,11)'ypos_verlet_noOO.dat'
  OPEN(unit=1,file=filename)
  WRITE(1,*) y
  CLOSE(1)


  DEALLOCATE(x,y,vx,vy,r)

END SUBROUTINE verlet
