PROGRAM test
  ! .exe compilation done by: $ gfortran -O3 -fopt-info -o task_b.exe task_b.f08 maxoffdiag_rotate.o
  USE maxoffdiag_rotate
  IMPLICIT NONE
  character(25) :: filename
  INTEGER :: i,j,k,l,n,iterations,u,p
  Real :: max,eps,h,c,s,testy,max_iterations
  double precision, allocatable, dimension(:) :: eig(:)
  double precision, allocatable, dimension(:,:) :: a(:,:),r(:,:)
  REAL, PARAMETER :: PI = 3.14159265358979323846264338327950288419

  ! Doloop variable
  p = 100

  ! Loop to calculate number of iterations for different matrix dimensions.
  do u = 1,p

    ! Variables
    n = u
    iterations = 0    ! Iteration counter
    eps = 10**(-8)    ! Accepted error
    h = 1.0/n     ! Step size

    ! Allocating
    allocate(a(n,n),r(n,n),eig(n))

    ! Setting up tridiagonal matrix
    do i = 1,n
      do j = 1,n
        if (i == j) then
          a(i,j) = 2.0/(h**2)
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

    ! First call to find max
    call maxoffdiag(n,a,max,l,k)

    ! Max number of iterations allowed.
    max_iterations = n*n*n

    ! Update(rotatetion) matrix a untill all non-diagonal elements are zero or
    ! iterations larger than max allowed iterations.
    do while (ABS(max) > eps .AND. iterations < max_iterations)
      call maxoffdiag(n,a,max,l,k)
      call rotate(n,a,r,l,k)
      iterations = iterations + 1
    enddo

    ! Writing to file can only be done after deleting existing dimension and iteration files as
    ! the appending just keeps writing!!
    10 format(A14)
    WRITE(filename,10)'iterations.dat'
    open(unit=1,file=filename,action='write',position='append')
    write(1,*) iterations
    close(1)

    11 format(A14)
    WRITE(filename,11)'dimensions.dat'
    open(unit=1,file=filename,action='write',position='append')
    write(1,*) n
    close(1)

    deallocate(a,r,eig)

enddo
ENDPROGRAM test
