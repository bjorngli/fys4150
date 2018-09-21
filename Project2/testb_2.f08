PROGRAM test
  IMPLICIT NONE
  character(25) :: filename
  INTEGER :: i,j,k,l,n,iterations,u,p
  Real :: max,eps,h,c,s,testy,max_iterations
  double precision, allocatable, dimension(:) :: eig(:)
  double precision, allocatable, dimension(:,:) :: a(:,:),r(:,:)
  REAL, PARAMETER :: PI = 3.14159265358979323846264338327950288419

  p = 100

  do u = 2,p
    ! Variables
    n = u
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

    ! Writing to file can only be done after deleting existing dimension and iteration files as
    ! the appending just keeps writing.
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




    !do i = 1,n
    !  eig(i) = 2.0/(h**2) + 2*(-1.0/(h**2))*cos(i*PI/(N+1.0))
    !enddo
    !write(*,*) eig
    deallocate(a,r,eig)

enddo
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
