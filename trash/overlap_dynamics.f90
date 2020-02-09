subroutine overlap_dynamics 
    use variables_definition 
    implicit none 
    real cap 
    integer :: i, k
    integer, parameter :: overlap_step = 1000
    real :: dt2, dt22, time
    
    dt2 = dt / 2.0 
    dt22 = dt * dt2 
    time = step * dt 

    cap = 0.001
    f(1:npart,1:dim) = cap + 0.05
    do k = 1, 1000
        call BOND_COND
        do i = 1, npart
            r(i,1:dim) = r(i,1:dim) + v(i,1:dim) * dt + dt22 * (f(i,1:dim)) / mass 
        enddo  
        !call FORCE 
        f(1:npart,1:dim) = f(1:npart,1:dim) + cap 
        do i = 1, npart 
           v(i,1:dim) =  v(i,1:dim) + dt2 * f(i,1:dim) / mass
        enddo 
        print*, f(1,:)
    end do 
    f(1:npart,1:dim) = 0.d0
    v(1:npart,1:dim) = 0.d0
    call ovito
    return 



end subroutine overlap_dynamics