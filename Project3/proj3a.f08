PROGRAM proj3a
IMPLICIT NONE
DOUBLE PRECISION :: FinalTime, h, pi, FourPi2
DOUBLE PRECISION, ALLOCATABLE :: time(:), x(:), y(:), vx(:), vy(:), r(:)
INTEGER :: n, i
CHARACTER(15) ::filename
FinalTime = 10000.0
n = 1000000
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

DO i=1,n+1
    time(i+1) = time(i) + h
    x(i+1) = x(i) + (h*vx(i))
    y(i+1) = y(i) + (h*vy(i))
    vx(i+1) = vx(i) - h*FourPi2*x(i+1)/(r(i)*r(i)*r(i))
    vy(i+1) = vy(i) - h*FourPi2*y(i+1)/(r(i)*r(i)*r(i))
    r(i+1) = sqrt((x(i+1)**2)+(y(i+1)**2))
    WRITE (*,*) time(i)
ENDDO

10 format(A15)
WRITE(filename,10)'xposition_a.dat'
open(unit=1,file=filename)
write(1,*) x
close(1)

11 format(A15)
WRITE(filename,11)'yposition_a.dat'
open(unit=1,file=filename)
write(1,*) y
close(1)

DEALLOCATE(time,x,y,vx,vy,r)
ENDPROGRAM
