PROGRAM test
  IMPLICIT NONE
  INTEGER :: i,j,k,l,n,iterations
  Real :: max,eps,h,c,s,testy,max_iterations
  double precision, allocatable, dimension(:) :: eig(:)
  double precision, allocatable, dimension(:,:) :: a(:,:),r(:,:)!,b(:,:)
  REAL, PARAMETER :: PI = 3.14159265358979323846264338327950288419
  real :: start, finish


    ! Variables
    n = 50
    iterations = 0
    eps = 10**(-8)
    h = 1./n

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

    call cpu_time(start)
    ! First call to find max
    call maxoffdiag(n,a,max,l,k)
    max_iterations = n*n*n

    ! Find max value not on diagonal and which array element,
    ! then do rotation which updates matrices a and r
    do while (ABS(max) > eps .AND. iterations < max_iterations)
      call maxoffdiag(n,a,max,l,k)
      call rotate(n,a,r,l,k)
      iterations = iterations + 1
    enddo
    call cpu_time(finish)
    !!WRITE(*,*) finish-start, 'sec'

    !write(*,*) a

    !do i = 1,n
    !  eig(i) = 2.0/(h**2) + 2*(-1.0/(h**2))*cos(i*PI/(N+1.0))
    !enddo
    !write(*,*) eig
    deallocate(a,r,eig)

    call unittest_maxoffdiag()
    call unittest_ortho_trans()

ENDPROGRAM test

subroutine maxoffdiag(num,mat,max,l,k)
  IMPLICIT NONE
  ! Input arguments
  integer, intent(in) :: num
  double precision, dimension(num,num), intent(in) :: mat
  ! Output arguments
  real, intent(out) :: max
  integer, intent(out) :: l,k
  ! Subroutine only arguments
  integer :: i,j
  max = 0.0

  do i = 1,num
   do j = i+1,num
     if (ABS(mat(i,j))>max) then
       max = ABS(mat(i,j))
       l = i
       k = j
     endif
   enddo
 enddo
end subroutine maxoffdiag

subroutine rotate(n,a,r,l,k)
  IMPLICIT NONE
  ! Input arguments
  integer, intent(in) :: n
  integer, intent(in) :: l,k
  ! In/out arguments
  double precision, dimension(n,n), intent(inout) :: a,r
  ! subroutine only arguments
  integer :: i
  real :: s,c,t,tau
  real :: a_kk,a_ll,a_ik,a_il,r_ik,r_il

  if (a(k,l) /= 0.0) then
    tau = (a(l,l) - a(k,k))/(2*a(k,l))
    if (tau > 0) then
      t = 1.0/(tau + sqrt(1.0+tau*tau))
    else
      t = -1.0/(-tau + sqrt(1.0+tau*tau))
    endif
    c = 1/sqrt(1+t*t)
    s = c*t
  else
    c = 1.0
    s = 0.0
  endif

! Take out values for later use (matrix gets rewritten so a_kk,a_ll would change)
a_kk = a(k,k)
a_ll = a(l,l)

a(k,k) = c*c*a_kk - 2.0*c*s*a(k,l) + s*s*a_ll
a(l,l) = s*s*a_kk + 2.0*c*s*a(k,l) + c*c*a_ll
a(k,l) = 0.0
a(l,k) = 0.0

do i = 1,n
  if (i /= k .AND. i /= l) then
    a_ik = a(i,k)
    a_il = a(i,l)
    a(i,k) = c*a_ik - s*a_il
    a(k,i) = a(i,k)
    a(i,l) = c*a_il + s*a_ik
    a(l,i) = a(i,l)
  endif
  r_ik = r(i,k)
  r_il = r(i,l)
  r(i,k) = c*r_ik - s*r_il
  r(i,l) = c*r_il + s*r_ik
enddo
end subroutine rotate

subroutine unittest_maxoffdiag()
! Unittest to check if algorithm searching for the largest non-diagonal element
! return the correct answer.

  IMPLICIT NONE
  ! In/out arguments
  double precision, allocatable, dimension(:,:) :: b(:,:)
  ! subroutine only arguments
  integer :: i,j,l,k,n
  real :: max

  n = 5

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


    ! Test 1. Check if a negative non-diagonal returns the correct value.
    ! Expect: max = 1 because max = ABS(-1).
    ! Check if difference between returned value and expected max equals 0
    ! within some margin.
    call maxoffdiag(n,b,max,l,k)
    if (ABS(max-(ABS(-1)))<=10**(-4)) then
      write(*,*) 'Max value test 1 successfull'!'Hello', max-ABS(-1), 10**(-5)
    else
      write(*,*) 'Max value test 1 failed'
    endif

    ! Test 2. Check if largest value on lower rows of matrix is returned when
    ! all values above is smaller.
    b(2,5) = -20
    !b(5,4) = -20
    !write(*,*) b
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
    if (ABS(max-(ABS(-20)))<=10**(-4)) then
      write(*,*) 'Max value test 2 successfull'!'Hello', max-ABS(-1), 10**(-5)
    else
      write(*,*) 'Max value test 2 failed'
    endif


end subroutine unittest_maxoffdiag

subroutine unittest_ortho_trans()
  ! Unit test to check if a matrix A has the same eigenvalues as after subroutine
  ! rotate is used to do orthogonal transformation.
  ! Done by comparing calculated eigenvalues of A by hand up against eigenvalues by
  ! lapack of matrix B where B is created by using a orthogonal transform on A

  ! Or using a simple matrix where you can calculate eigenvalues by hand.
  IMPLICIT NONE
  ! subroutine only arguments
  double precision, allocatable, dimension(:,:) :: a(:,:),r(:,:)
  integer :: i,j,n,l,k
  real :: max,eig_1,eig_2,eig_3

  n = 3

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
  !write(*,*) a

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

  ! Calculated eigenvalues by hand -> char. eq. r**2 - 4r + 2 = 0
  eig_1 = 2.0
  eig_2 = 2.0 - sqrt(2.0)
  eig_3 = 2.0 + sqrt(2.0)

  ! Calculated eigenvalues by jacobi's method rotation
  ! Eig values are on diagonal of matrix

  ! First call to find max
  call maxoffdiag(n,a,max,l,k)
  call rotate(n,a,r,l,k)
  write(*,*) a

end subroutine unittest_ortho_trans
