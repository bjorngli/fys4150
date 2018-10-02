PROGRAM project2b
IMPLICIT NONE
INTEGER :: i, j, N, LDZ, INFO
CHARACTER(1) :: JOBZ
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: Z(:,:), E(:), D(:), WORK(:)
DOUBLE PRECISION :: c, s, t, tau, theta, h
real :: start, finish
c = dcos(theta)
s = dsin(theta)
t = s/c
JOBZ = 'V'
WRITE(*,*) 'Enter a value for N: '
READ(*,*) N
LDZ = N
h = 1./N
ALLOCATE(E(N-1),Z(N,N),D(N),WORK(max(1,2*N-2)))

DO i=1,N
  E(i)=-1./(h**2)
  D(i)=2./(h**2)
ENDDO

call cpu_time(start)
CALL dstev(JOBZ,N,D,E,Z,LDZ,WORK,INFO)
call cpu_time(finish)

WRITE(*,*) finish-start, 'sec'
!WRITE(*,*) D, Z
DEALLOCATE(E,Z,D,WORK)
END PROGRAM
