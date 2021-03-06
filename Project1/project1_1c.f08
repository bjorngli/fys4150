program project1_1c
implicit NONE

character(22) :: filename
INTEGER :: num ! User input number
INTEGER :: i, n
REAL :: h
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: b(:),u(:),f(:),a(:),c(:) ! a and c arrays are the same
REAL :: start,finish

! Get n from user
  WRITE(*,*) 'Enter a number N: '
  READ(*,*) num
  n = 10**num
  ALLOCATE (a(n+1),b(n+1),f(n+1),u(n+1),c(n+1) ) ! Initializing vectors


  h = 1.0/n
  b(1)=2; b(n+1)=2; u(1)=0; u(n+1)=0
  ! Expression for f
  DO i=1,n+1
    f(i) = (h**2)*100*exp(-10*(i-1.0)*h) ! Many flops
  ENDDO

! Algorithm
  ! Updating b
  DO i=2,n
    b(i) = (i)/(i-1.0) ! 2 flops
  ENDDO

  call cpu_time(start)
  ! Forward substitution
  DO i=3,n
    f(i) = f(i) + f(i-1)/b(i-1) ! 2 flops
  ENDDO

  ! Backward substitution
  u(n) = f(n)/b(n)

  DO i=n-1,2,-1
    u(i) = (f(i)+u(i+1))/b(i) ! 2 flops
  ENDDO
! End Algorithm
  call cpu_time(finish)
  WRITE(*,*) finish-start, 'sec'

  10 format(A8,I1,A4)
  WRITE(filename,10)'s_vector',num,'.dat'

  open(unit=1,file=filename)
  write(1,*) u
  close(1)

  DEALLOCATE(a,b,c,f,u)
end program
