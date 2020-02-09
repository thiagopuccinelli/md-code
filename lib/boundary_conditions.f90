subroutine BOND_COND
    use variables_definition  
    integer :: i,k  

    do k = 1, dim 
        do i = 1, npart
            if (r(i,k) .lt. 0.d0) then
                r(i,k) = r(i,k) + BOX
            else if (r(i,k) .gt. BOX) then 
                r(i,k) = r(i,k) - BOX
            end if
        enddo 
    enddo 

    return 
end subroutine BOND_COND