program project1_1b
implicit none

character(21) :: filename
  INTEGER :: num ! User input number
  INTEGER :: i, n
  REAL :: h
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: b(:),u(:),f(:),a(:),c(:),frac(:) ! a and c arrays are the same
  REAL :: start,finish

  ! Get n from user
  WRITE(*,*) 'Enter a number N: '
  READ(*,*) num
  n = 10**num

    ALLOCATE(a(n+1),b(n+1),f(n+1),u(n+1),c(n+1),frac(n+1)) ! Initializing vectors


    h = 1.0/n
    b(1)=2; b(n+1)=2; u(1)=0; u(n+1)=0

    ! Equation for f
    DO i=1,n+1
      f(i) = (h**2)*100*exp(-10*(i-1.0)*h)
    ENDDO

    ! Setting up vectors
    DO i=1,n
      a(i) = -1
      c(i) = -1
      b(i) = 2
    ENDDO

    call cpu_time(start)
! Algorithm
    ! Update b
    DO i=3,n
      frac(i) = a(i-1)/b(i-1)      ! 1 flop
      b(i) = b(i)-(frac(i)*c(i-1)) ! 2 flops
      !Forward substitution
      f(i) = f(i) - f(i-1)*frac(i) ! 2 flops
    ENDDO

    ! Backward substitution
    u(n) = f(n)/b(n) ! 1 flop

    DO i=n-1,2,-1
      u(i) = (f(i)-c(i)*u(i+1))/b(i) ! 3 flops
    ENDDO
! End Algorithm
    call cpu_time(finish)
    WRITE(*,*) finish-start, 'sec'

    10 format(A8,I1,A4)
    WRITE(filename,10)'g_vector',num,'.dat'

    open(unit=1,file=filename)
    write(1,*) u
    close(1)

    DEALLOCATE(a,b,c,frac,f,u)
end program
