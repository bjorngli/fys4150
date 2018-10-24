module solver_class

  type solver
     integer::Bodies
     real(8),dimension(2) :: mass
     real(8),dimension(2,3) :: position
     real(8),dimension(2,3) :: velocity
  contains
    procedure ::relative_position
    procedure ::forces
    procedure ::calc_position
    procedure ::calc_velocities

 end type solver

contains
  subroutine relative_position(system,numbodies,relposition)
    implicit none

    class (solver),intent(in)::system
    integer :: i,j
    integer,intent(in)::NumBodies
    real(8),dimension(Numbodies,Numbodies,4):: relposition

    do i=1,NumBodies
       do j =1,NumBodies
          if (i.ne.j) then
             relposition(j,i,1:3)=system%position(j,1:3)-system%position(i,1:3)
             relposition(j,i,4)=sqrt(relposition(j,i,1)**2+relposition(j,i,2)**2+relposition(j,i,3)**2)
          end if
       end do
    end do

    return
  end subroutine relative_position

  subroutine forces(system,Numbodies,relposition,relforce,beta)
    implicit none

    class(solver),intent(in)::system
    real(8),intent(in)::beta
    integer :: j,i
    integer,intent(in):: Numbodies
    real(8),dimension(Numbodies,3)::relforce
    real(8),dimension(Numbodies,Numbodies,4):: relposition
    real(8)::rrr,Fourpi2,pow
    Fourpi2 = 4.d0*3.14*3.14
    pow = beta + 1

    do i=1,numbodies
       do j=1,3
          relforce(i,j)=0.d0
       end do
    end do


    do i=1,numbodies
       do j=1,numbodies
          if(j.ne.i) then
             rrr=(relposition(j,i,4)**pow)
             relforce(i,1:3) =relforce(i,1:3) - Fourpi2*system%mass(j)*relposition(j,i,1:3)/rrr
          end if
       end do
    end do
  end subroutine forces


  subroutine calc_position(system,Numbodies,relforce,h)
    implicit none
    class(solver), intent(inout)::system
    integer,intent(in)::Numbodies
    integer::i,j
    real(8),dimension(Numbodies,3)::relforce
    real(8)::h

    do i=1,numbodies
       system%position(i,1:3)=system%position(i,1:3)+h*system%velocity(i,1:3) - (h*h/2.d0)*relforce(i,:3)
    end do

  end subroutine calc_position

  subroutine calc_velocities(system,numbodies,relforce,updatedforce,h)
    implicit none
    class(solver),intent(inout)::system
    integer,intent(in)::NumBodies
    integer::i,j
    real(8)::h
    real(8),dimension(Numbodies,3)::relforce
    real(8),dimension(Numbodies,3)::updatedforce

    do i=1, numbodies
       system%velocity(i,1:3)=system%velocity(i,1:3)-h*0.5*updatedforce(i,1:3)-h*0.5*relforce(i,1:3)
    end do

  end subroutine calc_velocities

end module solver_class
