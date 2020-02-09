subroutine temperature_control
    use variables_definition  
    implicit none 
    integer :: i, k  
    real :: v2
    v2 = 0.d0 
    temperature = 0.d0
    en_kin = 0.d0
    do k = 1, dim  
        do i = 1, npart
            v2 = v2 + v(i,k) ** 2
        enddo 
    enddo 
    temperature  = v2 / (3*npart)
    v2 = v2 / npart 
    en_kin = 0.5d0 * mass * v2 
    open(30, file='ENERGY.dat')
    write(30, 3000) step, temperature, en_kin, sum(en_pot(:))/npart
3000 format(4X, I6, 3F20.5)
    return 
end subroutine temperature_control
