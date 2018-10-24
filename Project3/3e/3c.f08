program solarsystem
  use solver_class

  implicit none

  integer,parameter::numbodies=3
  real(8)::Tmax=50,h
  integer::Num_Steps=50000,m
  real(8),dimension(3)::masses
  real(8),dimension(3,3)::position
  real(8),dimension(3,3)::velocity
  real(8),dimension(numbodies,numbodies,4)::relposition,updaterel
  real(8),dimension(numbodies,3)::relforce
  real(8),dimension(numbodies,3)::updatedforce

  type(solver)::solar

  solar%mass(1)=1.
  solar%mass(2)=5.97219/1988500.
  solar%mass(3)=1898.13/1988500.

  solar%position(1,1)=-1.247339838016828E-04
  solar%position(1,2)=7.241110103144190E-03
  solar%position(1,3)=-7.324971134975035E-05
  solar%velocity(1,1)=365.25*(-7.578160673408990E-06)
  solar%velocity(1,2)=365.25*(2.617604518160534E-06)
  solar%velocity(1,3)=365.25*(1.894069341194622E-07)
  solar%position(2,1)=9.528047055398201E-01
  solar%position(2,2)=3.053612869840809E-01
  solar%position(2,3)=-9.272902073041313E-05
  solar%velocity(2,1)=365.25*(-5.428888690270241E-03)
  solar%velocity(2,2)=365.25*(1.636353485946535E-02)
  solar%velocity(2,3)=365.25*(-4.491683144318728E-07)
  solar%position(3,1)=-2.679859418467178E+00
  solar%position(3,2)=-4.648870533678862E+00
  solar%position(3,3)=7.922561878424454E-02
  solar%velocity(3,1)=365.25*(6.448353939132267E-03)
  solar%velocity(3,2)=365.25*(-3.409518873130117E-03)
  solar%velocity(3,3)=365.25*(-1.300875035507015E-04)

  open(2,file="3c_Sun.dat")
  open(3,file="3c_Earth.dat")
  open(4,file="3c_Jupiter.dat")

  h=Tmax/Num_Steps
  do m=1,Num_Steps

     call relative_position(solar,numbodies,relposition)
     call forces(solar,Numbodies,relposition,relforce)
     call calc_position(solar,numbodies,relforce,h)
     call relative_position(solar,numbodies,updaterel)
     call forces(solar,numbodies,updaterel,updatedforce)
     call calc_velocities(solar,numbodies,relforce,updatedforce,h)

     write(2,*)solar%position(1,1),solar%position(1,2),solar%position(1,3)
     write(3,*)solar%position(2,1),solar%position(2,2),solar%position(2,3)
     write(4,*)solar%position(3,1),solar%position(3,2),solar%position(3,3)
  end do

  close(2)
  close(3)
  close(4)

end program solarsystem
