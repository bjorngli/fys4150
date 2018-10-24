program solarsystem
  use solver_class

  implicit none

  integer,parameter::numbodies=2
  real(8)::Tmax=100,h
  integer::Num_Steps=100000,m
  real(8),dimension(2)::masses
  real(8),dimension(2,3)::position
  real(8),dimension(2,3)::velocity
  real(8),dimension(numbodies,numbodies,4)::relposition,updaterel
  real(8),dimension(numbodies,3)::relforce
  real(8),dimension(numbodies,3)::updatedforce

  ! Derived type object parameters.
  type(solver)::solar

  ! Initial conditions
  solar%mass(1)=1.d0
  solar%mass(2)=3.d-6

  solar%position(1,1)=0
  solar%position(1,2)=0
  solar%position(1,3)=0
  solar%velocity(1,1)=0
  solar%velocity(1,2)=0
  solar%velocity(1,3)=0
  solar%position(2,1)=1
  solar%position(2,2)=0
  solar%position(2,3)=0
  solar%velocity(2,1)=0
  solar%velocity(2,2)=2*3.14
  solar%velocity(2,3)=0

  ! Open file
  open(3,file="3b_Earth_verlet.dat")

  h=Tmax/Num_Steps
  do m=1,Num_Steps
     ! Call verlet method solver
     call relative_position(solar,numbodies,relposition)
     call forces(solar,Numbodies,relposition,relforce)
     call calc_position(solar,numbodies,relforce,h)
     call relative_position(solar,numbodies,updaterel)
     call forces(solar,numbodies,updaterel,updatedforce)
     call calc_velocities(solar,numbodies,relforce,updatedforce,h)
     ! Read to file
     write(3,*)solar%position(2,1),solar%position(2,2),solar%position(2,3)

  end do
  ! Close file
  close(3)

end program solarsystem
