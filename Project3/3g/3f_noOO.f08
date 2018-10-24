

Program verlet
  IMPLICIT NONE
  INTEGER :: n,i,outfile ! Number of planets
  real(8) :: solar_mass, earth_mass,FourPi2,FinalTime,h,solar_mass_value
  real(8) :: l, c=173*365.25,r_p,r,r1,x,x1,y,y1,vx,vy,a,a1
  CHARACTER(20) ::filename

  ! Variables
  FinalTime = 100.0
  n = 1e9
  h = FinalTime/REAL(n)
  FourPi2 = 4.d0*3.14*3.14
  solar_mass_value = 1.989e30
  solar_mass = 1 ! 1 AU
  earth_mass = 1.66012d-7 ! mercury mass in units of AU

  open(100,file="3f_pos.dat")

  ! Initial conditions
  x = 0.3075
  y = 0.0
  vx = 0.0
  vy = 12.44
  r = sqrt((x**2)+(y**2))
  r_p = x

  DO i=1,n
      l = abs(x*vy - y*vx)
      a = FourPi2/(r**3)*(1+(3*l**2)/(r**2*c**2))
      x1 = x + h*vx - ((h**2)/2)*(FourPi2)*x/(r**3)
      y1 = y + h*vy - ((h**2)/2)*(FourPi2)*y/(r**3)
      r1 = sqrt((x1**2)+(y1**2))
      a1 = FourPi2/(r1**3)*(1+(3*l**2)/(r1**2*c**2))
      vx = vx - (h/2)*(a1*x1+a*x)
      vy = vy - (h/2)*(a1*y1+a*y)
      if (r .lt. r1 .AND. r .lt. r_p) then
        WRITE(100,*)x,y
      end if
      r_p = r
      r = r1
      x = x1
      y = y1

  ENDDO

  CLOSE(100)

END PROGRAM verlet
