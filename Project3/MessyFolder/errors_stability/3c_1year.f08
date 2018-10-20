program test
  implicit None
  integer::i
  real(8)::dist_error,last_dist_error,good_n

  last_dist_error = 0
  good_n = 0
  open(300,file="dist_error.dat")
  do i=10,400
    call solarsystem(i,dist_error)
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
          write(*,*) good_n
        endif
      endif
    endif
  enddo
  close(300)

endprogram test

subroutine solarsystem(Num_Steps,dist_error)
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

end subroutine solarsystem
