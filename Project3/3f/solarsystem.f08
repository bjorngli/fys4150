program solarsystem
  use solver_class

  implicit none

  integer,parameter::numbodies=10
  real(8)::Tmax=250,h
  integer::Num_Steps=10000,m
  real(8),dimension(10)::masses
  real(8),dimension(10,3)::position
  real(8),dimension(10,3)::velocity
  real(8),dimension(numbodies,numbodies,4)::relposition,updaterel
  real(8),dimension(numbodies,3)::relforce
  real(8),dimension(numbodies,3)::updatedforce
  !real(8),dimension(numbodies+1)::kinetic,potential,angular

  type(solver)::solar


  solar%mass(1)=1.                !Sun
  solar%mass(2)=5.97219/1988500.  !Earth
  solar%mass(3)=1898.13/1988500.  !Jupiter
  solar%mass(4)=0.64171/1988500.  !Mars
  solar%mass(5)=4.8685/1988500.   !Venus
  solar%mass(6)=568.34/1988500.   !Saturn
  solar%mass(7)=0.3302/1988500.   !Mercury
  solar%mass(8)=86.813/1988500.   !Uranus
  solar%mass(9)=102.413/1988500.  !Neptune
  solar%mass(10)=0.01307/1988500. !Pluto

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
  solar%position(4,1)=1.370566065802587E+00
  solar%position(4,2)=-1.928843729924704E-01
  solar%position(4,3)=-3.790047470953376E-02
  solar%velocity(4,1)=365.25*(2.548393097798315E-03)
  solar%velocity(4,2)=365.25*(1.504516606737405E-02)
  solar%velocity(4,3)=365.25*(2.526759216773895E-04)
  solar%position(5,1)=7.177791501659343E-01
  solar%position(5,2)=1.078222036794456E-01
  solar%position(5,3)=-4.012084731109569E-02
  solar%velocity(5,1)=365.25*(-2.892213210424378E-03)
  solar%velocity(5,2)=365.25*(1.994145100143018E-02)
  solar%velocity(5,3)=365.25*(4.401915553550370E-04)
  solar%position(6,1)=1.533542535074663E+00
  solar%position(6,2)=-9.937711138425550E+00
  solar%position(6,3)=1.117470354327505E-01
  solar%velocity(6,1)=365.25*(5.206732001657008E-03)
  solar%velocity(6,2)=365.25*(8.336898877431193E-04)
  solar%velocity(6,3)=365.25*(-2.220848255804162E-04)
  solar%position(7,1)=-2.138801978256535E-01
  solar%position(7,2)=-4.028538544608614E-01
  solar%position(7,3)=-1.397398026440086E-02
  solar%velocity(7,1)=365.25*(1.927259979627735E-02)
  solar%velocity(7,2)=365.25*(-1.164174161133437E-02)
  solar%velocity(7,3)=365.25*(-2.720041450765680E-03)
  solar%position(8,1)=1.718193423311543E+01
  solar%position(8,2)=9.988014846393089E+00
  solar%position(8,3)=-1.854987455727011E-01
  solar%velocity(8,1)=365.25*(-2.005621300007712E-03)
  solar%velocity(8,2)=365.25*(3.217042846869565E-03)
  solar%velocity(8,3)=365.25*(3.810139183912805E-05)
  solar%position(9,1)=2.891793183858609E+01
  solar%position(9,2)=-7.731693744208183E+00
  solar%position(9,3)=-5.072234730681663E-01
  solar%velocity(9,1)=365.25*(7.899417770798214E-04)
  solar%velocity(9,2)=365.25*(3.050877562321204E-03)
  solar%velocity(9,3)=365.25*(-8.139357296428597E-05)
  solar%position(10,1)=1.163483999459982E+01
  solar%position(10,2)=-3.157641175938865E+01
  solar%position(10,3)=1.338119585887108E-02
  solar%velocity(10,1)=365.25*(3.011486331290935E-03)
  solar%velocity(10,2)=365.25*(4.287174206963734E-04)
  solar%velocity(10,3)=365.25*(-9.044267681804708E-04)

  open(3,file="Earth.dat")
  open(4,file="Jupiter.dat")
  !open(5,file="energy.dat")
  !open(7,file="momentum.dat")
  open(8,file="Mars.dat")
  open(9,file="Venus.dat")
  open(10,file="Saturn.dat")
  open(11,file="Mercury.dat")
  open(12,file="Uranus.dat")
  open(13,file="Neptune.dat")
  open(14,file="Pluto.dat")


  h=Tmax/Num_Steps
  do m=1,Num_Steps

     call relative_position(solar,numbodies,relposition)
     call forces(solar,Numbodies,relposition,relforce)
     call calc_position(solar,numbodies,relforce,h)
     call relative_position(solar,numbodies,updaterel)
     call forces(solar,numbodies,updaterel,updatedforce)
     call calc_velocities(solar,numbodies,relforce,updatedforce,h)
     !call kinetic_energy(solar,numbodies,kinetic)
     !call potential_energy(solar,numbodies,updaterel,potential)
     !call angular_momentum(solar,numbodies,updaterel,angular)

!     write(5,*), kinetic(Numbodies+1), potential(numbodies+1)
!     write(7,*), angular(Numbodies+1)

     write(3,*)solar%position(2,1),solar%position(2,2),solar%position(2,3)
     write(4,*)solar%position(3,1),solar%position(3,2),solar%position(3,3)
     write(8,*)solar%position(4,1),solar%position(4,2),solar%position(4,3)
     write(9,*)solar%position(5,1),solar%position(5,2),solar%position(5,3)
     write(10,*)solar%position(6,1),solar%position(6,2),solar%position(6,3)
     write(11,*)solar%position(7,1),solar%position(7,2),solar%position(7,3)
     write(12,*)solar%position(8,1),solar%position(8,2),solar%position(8,3)
     write(13,*)solar%position(9,1),solar%position(9,2),solar%position(9,3)
     write(14,*)solar%position(10,1),solar%position(10,2),solar%position(10,3)

  end do


  close(3)
  close(4)
  !close(5)
  !close(7)
  close(8)
  close(9)
  close(10)
  close(11)
  close(12)
  close(13)
  close(14)
end program solarsystem
