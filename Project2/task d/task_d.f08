PROGRAM test
  ! .exe compilation done by: $ gfortran -O3 -fopt-info -o test_d.exe test_d.f08 maxoffdiag_rotate.o
  USE maxoffdiag_rotate
  IMPLICIT NONE
  INTEGER :: i,j,k,l,n,iterations
  Real :: max,eps,h,c,s,testy,max_iterations,hh,rho_N
  double precision, allocatable, dimension(:) :: rho(:),eig_one
  double precision, allocatable, dimension(:,:) :: a(:,:),r(:,:)
  REAL, PARAMETER :: PI = 3.14159265358979323846264338327950288419
  real :: start, finish
  character(25) :: filename


    ! Variables
    n = 200
    iterations = 0
    eps = 10**(-8)
    rho_N = 5
    h = rho_N/n

    ! Allocating
    allocate(a(n,n),r(n,n),rho(n),eig_one(n))

    ! Potential array
    do i = 1,n
      rho(i) = i*h
    enddo


    ! Setting up tridiagonal matrix
    do i = 1,n
      do j = 1,n
        if (i == j) then
          a(i,j) = (2.0/(h**2)) + rho(i)**2
        elseif (i==(j+1)) then
          a(i,j) = -1.0/(h**2)
        elseif (i==(j-1)) then
          a(i,j) = -1.0/(h**2)
        else
          a(i,j) = 0
        endif
      end do
    end do

    ! Setting up eigenvector matrix
    do i = 1,n
      do j = 1,n
          if (i == j) then
            r(i,j) = 1.0
          else
            r(i,j) = 0.0
          endif
      enddo
    enddo

    ! Timecheck
    call cpu_time(start)

    ! First call to find max
    call maxoffdiag(n,a,max,l,k)
    max_iterations = n*n*n

    ! Find max value not on diagonal and which array element
    ! then do rotation which updates matrices a and r
    do while (ABS(max) > eps .AND. iterations < max_iterations)
      call maxoffdiag(n,a,max,l,k)
      call rotate(n,a,r,l,k)
      iterations = iterations + 1
    enddo

    ! End timecheck and write
    call cpu_time(finish)
    WRITE(*,*) finish-start, 'sec'

    ! Send eigenvalues to array eig_one and write to file
    do i = 1,n
      do j = 1,n
        if (i==j) then
          eig_one(i) = a(i,j)
        endif
      enddo
    enddo

    deallocate(a)

    10 format(A16)
    WRITE(filename,10)'one_electron.dat'
    open(unit=1,file=filename)
    write(1,*) eig_one
    close(1)

    deallocate(eig_one)

    open (unit=2, file='one_outfile.dat')
    do  i=1,n
        write(2,*) r(i,:)
    end do
    close(2)

    deallocate(r)

ENDPROGRAM test
