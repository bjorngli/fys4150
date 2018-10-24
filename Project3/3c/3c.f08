! Run example
! $ gfortran -O3 -fopt-info -o 3c 3c.f08 3c_solver.o

program pos
  call error_euler()
  call error_verlet()
  call energy_mom()

end program pos

subroutine error_euler()
  IMPLICIT None
  integer :: i,good_n
  DOUBLE PRECISION:: dist_error,last_dist_error

  last_dist_error = 0
  good_n = 0

  open(300,file="dist_error_euler.dat")
  do i=20,1500
    CALL euler(i,dist_error)
    write(300,*) dist_error
    if (good_n > 0) then
      continue
    else
      ! For first loop run, store current_dist_error
      if (i == 10) then
        last_dist_error = dist_error
      ! For all other loop runs
      else
        ! If change in distance error larger than 0.1%, keep searching
        if (abs(dist_error-last_dist_error)/last_dist_error*100 > 0.1) then
          last_dist_error = dist_error
        ! If error below 0.1%, store loop run as good enough amount of integration points.
        else
          good_n = i
          write(*,*) 'Good number of integration points:',good_n
        endif
      endif
    endif
  enddo
  close(300)

END subroutine error_euler

SUBROUTINE euler(n,dist_error)
  IMPLICIT NONE
  DOUBLE PRECISION :: FinalTime, h, pi, FourPi2,x_diff,y_diff
  DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:):: time, x, y, vx, vy
  DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:):: r
  INTEGER,intent(in) :: n
  DOUBLE PRECISION,intent(out)::dist_error
  INTEGER::i
  CHARACTER(19) ::filename
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: posx,posy,velx,vely,ax,ay,Fx,Fy

  FinalTime = 1.0
  h = FinalTime/REAL(n)
  pi = 4.D0*DATAN(1.D0)
  FourPi2 = 4*pi*pi
  ALLOCATE(time(n+1),x(n+1),y(n+1),vx(n+1),vy(n+1),r(n+1))
  time(1) = 0.0
  x(1) = 1.0
  y(1) = 0.0
  vx(1) = 0.0
  vy(1) = 2.0*pi
  r(1) = sqrt((x(1)**2)+(y(1)**2))

  DO i=1,n
      time(i+1) = time(i) + h
      x(i+1) = x(i) + (h*vx(i))
      y(i+1) = y(i) + (h*vy(i))
      vx(i+1) = vx(i) - h*FourPi2*x(i+1)/(r(i)*r(i)*r(i))
      vy(i+1) = vy(i) - h*FourPi2*y(i+1)/(r(i)*r(i)*r(i))
      r(i+1) = sqrt((x(i+1)**2)+(y(i+1)**2))
  ENDDO
  x_diff = x(1)-x(n)
  y_diff = y(1)-y(n)
  dist_error = sqrt((x_diff**2)+(y_diff**2))
  DEALLOCATE(time,x,y,vx,vy,r)
END SUBROUTINE euler

subroutine error_verlet()
  implicit None
  integer::i
  real(8)::dist_error,last_dist_error,good_n

  last_dist_error = 0
  good_n = 0
  open(300,file="dist_error.dat")
  do i=20,1500
    call verlet(i,dist_error)
    write(300,*) dist_error
    if (good_n > 0) then
      continue
    else
      ! For first loop run, store current_dist_error
      if (i == 10) then
        last_dist_error = dist_error
      ! For all other loop runs
      else
        ! If change in distance error larger than 0.1%, keep searching
        if (abs(dist_error-last_dist_error)/last_dist_error*100 > 0.1) then
          last_dist_error = dist_error
        ! If error below 0.1%, store loop run as good enough amount of integration points.
        else
          good_n = i
          write(*,*) 'Good number of integration points:', good_n
        endif
      endif
    endif
  enddo
  close(300)

end subroutine error_verlet

subroutine verlet(Num_Steps,dist_error)
  use solver_class

  implicit none

  integer, intent(in)::Num_Steps
  integer,parameter::numbodies=2
  real(8)::Tmax=1,h
  integer::m
  real(8),dimension(3)::masses
  real(8),dimension(3,3)::position
  real(8),dimension(3,3)::velocity
  real(8),dimension(numbodies,numbodies,4)::relposition,updaterel
  real(8),dimension(numbodies,3)::relforce
  real(8),dimension(numbodies,3)::updatedforce
  real(8),dimension(3)::pos_start
  real(8),dimension(3)::dist_error_vec
  real(8),intent(out)::dist_error
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

  h=Tmax/Num_Steps
  do m=1,Num_Steps
     if (m==1) then
       pos_start = solar%position(2,1:3)
     endif
     call relative_position(solar,numbodies,relposition)
     call forces(solar,Numbodies,relposition,relforce)
     call calc_position(solar,numbodies,relforce,h)
     call relative_position(solar,numbodies,updaterel)
     call forces(solar,numbodies,updaterel,updatedforce)
     call calc_velocities(solar,numbodies,relforce,updatedforce,h)
     if (m==Num_Steps) then
       dist_error_vec = pos_start - solar%position(2,1:3)
       dist_error = sqrt(dist_error_vec(1)**2+dist_error_vec(2)**2+dist_error_vec(3)**2)
     endif
  end do

end subroutine verlet

subroutine energy_mom()
  implicit None
  integer::i
  real(8)::energy_error,momentum_error

  open(300,file="energy_error.dat")
  open(301,file="momentum_error.dat")

  do i=20,500
    call energy_momentum(i,energy_error,momentum_error)

    ! After run of i timesteps write largest recorded value to file.
    write(300,*) energy_error
    write(301,*) momentum_error
  enddo
  close(300)
  close(301)
end subroutine energy_mom

subroutine energy_momentum(i,energy_error,momentum_error)
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
         ! If energy difference outside limit,
           ! If energy value stored is larger than current error. Skip. Else, record new max energy error
         if (abs(tot_energy-(kinetic(numbodies+1) + potential(numbodies+1))) >= 1e-8) then
           if (energy_error > abs(tot_energy-(kinetic(numbodies+1) + potential(numbodies+1)))) then
             continue
           else
             energy_error = abs(tot_energy-(kinetic(numbodies+1) + potential(numbodies+1)))
           endif
         endif
         ! If mom. difference outside limit
           ! If mom. value stored is larger than current error. Skip. Else, record new max mom. error
         if (abs((tot_ang)-(angular(numbodies+1))) >= 1e-10) then
           if (momentum_error > abs((tot_ang)-(angular(numbodies+1)))) then
             continue
           else
             momentum_error = abs((tot_ang)-(angular(numbodies+1)))
           endif
         endif
       endif
    end do

end subroutine energy_momentum
