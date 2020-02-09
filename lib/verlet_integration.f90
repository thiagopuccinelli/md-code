subroutine MOVE_VERLET
! 	** VELOCITY VERLET ALGORITHM **
!	** AT THE START OF A TIMESTEP, MOVE_VELVERLET_A IS CALLED TO ADVANCE THE 
!	** THE POSITIONS AND 'HALF-ADVANCE' THE VELOCITIES. THEN THE FORCE ROUTINE 
!	** IS CALLED, AND THIS IS FOLLOWED BY MOVE_VELVERLET_B WHICH COMPLETES THE 
!	** ADVANCEMENT OF VELOCITIES 
    use variables_definition 
    implicit none 
    integer :: i, k
    real :: dt2, dt22, time
    
    dt2 = dt / 2.0 
    dt22 = dt * dt2 
    time = step * dt 

    do i = 1, npart
        r(i,1:dim) = r(i,1:dim) + v(i,1:dim) * dt + dt22 * (f(i,1:dim)) / mass 
    enddo  
    call FORCE
    do i = 1, npart 
        v(i,1:dim) =  v(i,1:dim) + dt2 * f(i,1:dim) / mass
    enddo 
    
    return 
end subroutine MOVE_VERLET