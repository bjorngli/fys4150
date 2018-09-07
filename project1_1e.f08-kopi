program project1_1e
use F90library
implicit NONE
character(22) :: filename
double precision, allocatable, dimension(:,:) :: a(:,:)
double precision, allocatable, dimension(:) ::  f(:),v(:)
integer :: n, i, j, num
integer, allocatable :: indx(:)
real :: h, start, finish
double precision :: d

WRITE(*,*) 'Enter a value for N: '
READ(*,*) num
n = 10**num
allocate (a(n-1,n-1),f(n-1),indx(n-1),v(n+1))

h = 1.0/(n+1)

v(1)=0
v(n+1)=0
DO i=1,n-1
  f(i) = (h**2)*100*exp(-10*(i)*h)
ENDDO


do i = 1,n-1
  do j = 1,n-1
    if (i == j) then
      a(i,j) = 2
    elseif (i==(j+1)) then
      a(i,j) = -1
    elseif (i==(j-1)) then
      a(i,j) = -1
    else
      a(i,j) = 0
    endif
  end do
end do
call cpu_time(start)
call lu_decompose(a,n-1,indx,d)
call lu_linear_equation(a,n-1,indx,f)
call cpu_time(finish)

WRITE(*,*) finish-start, 'sec'
  do i = 1,n
    v(i+1)=f(i)
  end do

10 format(A9,I1,A4)
WRITE(filename,10)'lu_vector',num,'.dat'

open(unit=1,file=filename)
write(1,*) v
close(1)

end program
