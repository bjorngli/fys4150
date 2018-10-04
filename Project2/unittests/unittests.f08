PROGRAM unittests
  ! .exe compilation done by: $ gfortran -O3 -fopt-info -o unittests.exe unittests.f08 maxoffdiag_rotate.o
  IMPLICIT NONE
  call unittest_maxoffdiag()
  call unittest_orthogonal_transformation()
END PROGRAM unittests

subroutine unittest_maxoffdiag()
  ! Unittest to check if algorithm searching for the largest non-diagonal element
  ! return the correct answer.
  USE maxoffdiag_rotate
  IMPLICIT NONE
  ! In/out arguments
  double precision, allocatable, dimension(:,:) :: b(:,:)
  ! Subroutine only arguments
  integer :: i,j,l,k,n,tests
  real :: max,larg_val
  double precision :: epsillon

  ! Variables
  n = 5               ! Matrix dimension
  larg_val = -20.0    ! New largest value.
  tests = 0           ! Number of tests should equal 2 after code if 2 tests pass.
  epsillon = 10e-8    ! Accepted error.

  ! Allocate
  allocate(b(n,n))

  ! Setting up tridiagonal matrix
  do i = 1,n
    do j = 1,n
      if (i == j) then
        b(i,j) = 2
      elseif (i==(j+1)) then
        b(i,j) = -1.0
      elseif (i==(j-1)) then
        b(i,j) = -1.0
      else
        b(i,j) = 0
      endif
    end do
  enddo


    ! Test negative value
    ! Check if a negative non-diagonal returns the correct value.
    ! Expect: max = 1 because max = ABS(-1).
    ! Check if difference between returned value and expected max equals 0
    ! within some margin.
    call maxoffdiag(n,b,max,l,k)
    if (ABS(max-(ABS(-1)))<=epsillon) then
      !write(*,*) 'Max value test 1 successfull'!'Hello', max-ABS(-1), 10**(-5)
      tests = tests + 1
    else
      write(*,*) 'Subroutine maxoffdiag: Test 1 failed'
    endif

    ! Test new largest value
    ! Check if largest value on lower rows of matrix is returned when
    ! all values above is smaller.
    b(2,5) = larg_val

    ! maxoffdiag only looks at lower diagonal. Is this okay or should it look on upper
    ! digonal aswell?
    ! Also can one assume same values on lower and upper diogonal?
    ! maxoffdiag partial code:
    ! do i = 1,n
    !  do j = i+1,n
    !    write(*,*) b(i,j)
       ! if (ABS(b(i,j))>max) then
       !   max = ABS(b(i,j))
    !  enddo
    ! enddo
    call maxoffdiag(n,b,max,l,k)
    if (ABS(max-(ABS(-20)))<=epsillon) then
      !write(*,*) 'Max value test 2 successfull'!'Hello', max-ABS(-1), 10**(-5)
      tests = tests + 1
    else
      write(*,*) 'Subroutine maxoffdiag: Test 2 failed'
    endif

    if (tests == 2) then
      write(*,*) 'Subroutine maxoffdiag: All tests passed'
    endif
end subroutine unittest_maxoffdiag

subroutine unittest_orthogonal_transformation()
  ! Unit test to check if:
  ! 1. Eigenvectors are orthogonal before and after rotation (orthogonal transformation)
  ! 2. Matrix A has the same eigenvalues before and after rotation (orthogonal transformation)

  ! 1. Done by knowing the dot product of orthogonal vectors should equal zero
  ! 2. Done by comparing calculated eigenvalues of A by hand up against eigenvalues by
  ! lapack of matrix B where B is created by using a orthogonal transform on A
  ! Comment to self: Should maybe use hand calculation on both or lapack on both?
  USE maxoffdiag_rotate
  IMPLICIT NONE
  ! Subroutine only arguments
  double precision, allocatable, dimension(:,:) :: a(:,:),r(:,:)!,lapack
  integer :: i,j,n,l,k,tests
  real :: max,eig_1,eig_2,eig_3,lapack_1,lapack_2,lapack_3

  ! Variables
  n = 3         ! Matrix dimension
  tests = 0     ! Number of tests should equal 2 after code if 2 test pass

  ! Allocate
  allocate(a(n,n),r(n,n))

  ! Setting up tridiagonal matrix
  do i = 1,n
    do j = 1,n
      if (i == j) then
        a(i,j) = 2
      elseif (i==(j+1)) then
        a(i,j) = -1.0
      elseif (i==(j-1)) then
        a(i,j) = -1.0
      else
        a(i,j) = 0.0
      endif
    end do
  enddo

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


  ! Test orthogonality prework
  ! Make sure initial eigenvector is orthogonal. If not orthogonal write error message.
  if (DOT_PRODUCT(r(:,1),r(:,2)) + DOT_PRODUCT(r(:,1),r(:,3)) + DOT_PRODUCT(r(:,2),r(:,3))<=10e-5) then
      ! Test preserved eigenvalues prework
      ! Calculated eigenvalues by hand -> char. eq. r**2 - 4r + 2 = 0
      eig_1 = 2.0
      eig_2 = 2.0 - sqrt(2.0)
      eig_3 = 2.0 + sqrt(2.0)

      ! First call to find max and then rotate
      call maxoffdiag(n,a,max,l,k)
      ! Doing one orthogonal transformation
      call rotate(n,a,r,l,k)

      ! Test orthogonal eigenvectors
      if (DOT_PRODUCT(r(:,1),r(:,2)) + DOT_PRODUCT(r(:,1),r(:,3)) + DOT_PRODUCT(r(:,2),r(:,3))<=10e-5) then
        tests = tests + 1
      else
        write(*,*) 'Subroutine rotate: Orthogonality test failed'
      endif

      ! Test eigenvalues
      ! Insert LAPACK solving of new matrix A.
      ! Decide if you wanna do lapack as a matrix and index lapack(1) or do lapack_1,lapack_2... as variables.
      ! Just remember to put the variables/arrays under IMPLICIT NONE.
      !  if ( (eig_1 - lapack_1) <= 10**(-8) .AND. (eig_2 - lapack_2) <= 10**(-8) .AND. (eig_3 - lapack_3) <= 10**(-8)) then
      if ( (1.0 - 1.0) <= 10**(-8) .AND. (1.0 - 1.0) <= 10**(-8) .AND. (1.0 - 1.0) <= 10**(-8)) then
        !write(*,*) 'test'
        tests = tests + 1
      else
        write(*,*) 'Subroutine rotate: Eigenvalue test failed'
      endif

  else
    write(*,*) 'ERROR: Your initial eigenvectors are not orthogonal'
  endif

  if (tests == 2) then
    write(*,*) 'Subroutine rotate: All tests passed'
  endif
end subroutine unittest_orthogonal_transformation
