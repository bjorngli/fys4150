program solarsystem
  use solver_class

  implicit none

  integer,parameter::numbodies=2
  real(8)::Tmax=100,h
  integer::Num_Steps=100000,m,k
  real(8),dimension(3)::masses
  real(8),dimension(3,3)::position
  real(8),dimension(3,3)::velocity
  real(8),dimension(numbodies,numbodies,4)::relposition,updaterel
  real(8),dimension(numbodies,3)::relforce
  real(8),dimension(numbodies,3)::updatedforce
  real(8),dimension(5)::beta=[2.0,2.25,2.50,2.75,3.0]
  CHARACTER(20)::filename

  type(solver)::solar


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
  solar%velocity(2,2)=2*3.14   ! 2.8*3.14 = escape; 1*3.14 = elipse
  solar%velocity(2,3)=0
  do k =1,5
    10 FORMAT(A6,F4.2,A4)
    WRITE(filename,10)'earth_',beta(k),".dat"
    open(3+k,file=filename)

    h=Tmax/Num_Steps
    do m=1,Num_Steps

       call relative_position(solar,numbodies,relposition)
       call forces(solar,Numbodies,relposition,relforce,beta(k))
       call calc_position(solar,numbodies,relforce,h)
       call relative_position(solar,numbodies,updaterel)
       call forces(solar,numbodies,updaterel,updatedforce,beta(k))
       call calc_velocities(solar,numbodies,relforce,updatedforce,h)
       write(3+k,*)solar%position(2,1),solar%position(2,2),solar%position(2,3)

    end do
    close(3+k)

  end do

end program solarsystem
