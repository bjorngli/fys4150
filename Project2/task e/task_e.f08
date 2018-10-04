PROGRAM two_electron
  ! .exe compilation done by: $ gfortran -O3 -fopt-info -o test_e.exe test_e.f08 maxoffdiag_rotate.o

  ! Solving the scaled radial Schroedinger’s equation containing
  ! omega_r/freq as a parameter reflecting the strenght of the oscillator
  ! potential for four different values of omega_r/freq.
  ! For the ground state l = 0.
  ! Program solves by doing orthogonal transforms on a matrix created by
  ! taylor expanding Schroedinger’s equation and finding its relation for
  ! setting the equation up as a matrix problem Ax = λx, where A is a
  ! tridiagonal matrix.
  ! Then do orthogonal transforms until the non-diagonal elements are zero.
  ! What remains in matrix A are then the eigenvalues and r contains the
  ! eigenvectors.
  USE maxoffdiag_rotate
  IMPLICIT NONE
  INTEGER :: i,j,k,l,n,iterations,e
  Real :: max,eps,h,c,s,testy,max_iterations,hh,rho_N
  double precision, allocatable, dimension(:) :: rho(:),eig_one(:),freq(:)
  double precision, allocatable, dimension(:,:) :: a(:,:),r(:,:)
  REAL, PARAMETER :: PI = 3.14159265358979323846264338327950288419
  real :: start, finish
  character(25) :: filename

! Do loops runs code 4 times with different freq/omega_r each time.
do e = 1,4
    ! Variables
    n = 200             ! Matrix dimension
    iterations = 0      ! Iteration counter
    eps = 10e-8         ! Epsillon for checking error is not too big
    rho_N = 5           ! rho_max value
    h = rho_N/n         ! Calculting

    ! Allocate a,n,rho,eig_one
    allocate(a(n,n),r(n,n),rho(n),eig_one(n))
    ! Allocate four frequencies.
    allocate(freq(4))

    ! Sett up rho
    do i = 1,n
      rho(i) = i*h
    enddo

    ! Sett frequency values
    freq(1) = 0.01
    freq(2) = 0.5
    freq(3) = 1.0
    freq(4) = 5.0

    ! Sett up tridiagonal matrix
    do i = 1,n
      do j = 1,n
        if (i == j) then
          a(i,j) = (2.0/(h**2)) + (freq(e)**2)*(rho(i)**2) + 1.0/rho(i)
        elseif (i==(j+1)) then
          a(i,j) = -1.0/(h**2)
        elseif (i==(j-1)) then
          a(i,j) = -1.0/(h**2)
        else
          a(i,j) = 0
        endif
      end do
    end do

    ! Sett up eigenvector matrix
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
    ! Sett max iterations as n**3
    max_iterations = n*n*n

    ! Do while max non-diagonal value not close to zero. Or iter < n**3
    ! Find max value not on diagonal and which array element.
    ! Then, do rotation which updates matrices a and r.
    do while (ABS(max) > eps .AND. iterations < max_iterations)
      call maxoffdiag(n,a,max,l,k)
      call rotate(n,a,r,l,k)
      iterations = iterations + 1
    enddo

    ! Put eigenvalues in array eig_one and write to file.
    do i = 1,n
      do j = 1,n
        if (i==j) then
          eig_one(i) = a(i,j)
        endif
      enddo
    enddo

    10 format(A12,I1,A4)
    WRITE(filename,10)'two_electron',e,'.dat'
    open(unit=1,file=filename)
    write(1,*) eig_one
    close(1)

    11 format(A11,I1,A4)
    WRITE(filename,11)'two_outfile',e,'.txt'
    open (unit=2, file=filename)
    do  i=1,n
        write(2,*) r(i,:)
    end do

    ! Deallocate
    deallocate(a,r,eig_one,rho,freq)

enddo
ENDPROGRAM two_electron
