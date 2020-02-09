subroutine system_initialization 
    use variables_definition 
    implicit none 
    integer :: i, k 
    integer :: iseed
    real vel
    real fs, v2 !for the setting zero center of mass
    real, parameter :: temp = 1.d0 
    real gasdev 
    iseed = 387489
    r = 0 
    f = 0 
    v = 0 
    en_pot = 0.d0 
    en_kin = 0.d0
    call LATTICE
    do i = 1, npart 
        do k = 1, dim 
            vel = gasdev(iseed)
            v(i,k) = vel
        enddo 
    enddo 
    return 
end subroutine system_initialization 