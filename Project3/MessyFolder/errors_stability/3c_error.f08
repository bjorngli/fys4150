program test
  implicit None
  integer::i
  real(8)::energy_error,momentum_error

  open(300,file="energy_error.dat")
  open(301,file="momentum_error.dat")

  do i=20,500
    call solarsystem(i,energy_error,momentum_error)

    ! After run of i timesteps write largest recorded value to file.
    write(300,*) energy_error
    write(301,*) momentum_error
  enddo
  close(300)
  close(301)
endprogram test

subroutine solarsystem(i,energy_error,momentum_error)
  use solver_class

  implicit none

  integer, intent(in)::i
  integer,parameter::numbodies=2
  real(8)::Tmax=10,h
  integer::m
  real(8),dimension(3)::masses
  real(8),dimension(3,3)::position
  real(8),dimension(3,3)::velocity
  real(8),dimension(numbodies,numbodies,4)::relposition,updaterel
  real(8),dimension(numbodies,3)::relforce
  real(8),dimension(numbodies,3)::updatedforce
  real(8),dimension(numbodies+1)::kinetic,potential,angular
  real(8) :: tot_energy,tot_ang
  real(8),intent(out)::energy_error,momentum_error

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

    energy_error = 0
    momentum_error = 0
    h=Tmax/i
    do m=1,i

       call relative_position(solar,numbodies,relposition)
       call forces(solar,numbodies,relposition,relforce)
       call calc_position(solar,numbodies,relforce,h)
       call relative_position(solar,numbodies,updaterel)
       call forces(solar,numbodies,updaterel,updatedforce)
       call calc_velocities(solar,numbodies,relforce,updatedforce,h)
       call kinetic_energy(solar,numbodies,kinetic)
       call potential_energy(solar,numbodies,updaterel,potential)
       call angular_momentum(solar,numbodies,updaterel,angular)

       ! Store first energy and momentum record to compare against
       if (m==1) then
         tot_energy = kinetic(numbodies+1) + potential(numbodies+1)
         tot_ang = angular(numbodies+1)
       ! For all other m:
       else
         ! If energy difference outside limit, if energy value stored is larger than current error. Skip. Else, record new max energy error
         if (abs(tot_energy-(kinetic(numbodies+1) + potential(numbodies+1))) >= 1e-8) then
           if (energy_error > abs(tot_energy-(kinetic(numbodies+1) + potential(numbodies+1)))) then
             continue
           else
             energy_error = abs(tot_energy-(kinetic(numbodies+1) + potential(numbodies+1)))
           endif
         endif
         ! If mom. difference outside limit, if mom. value stored is larger than current error. Skip. Else, record new max mom. error
         if (abs((tot_ang)-(angular(numbodies+1))) >= 1e-10) then
           if (momentum_error > abs((tot_ang)-(angular(numbodies+1)))) then
             continue
           else
             momentum_error = abs((tot_ang)-(angular(numbodies+1)))
           endif
         endif
       endif
    end do

end subroutine solarsystem
