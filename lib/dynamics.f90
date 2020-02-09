subroutine dynamics 
    use variables_definition 
    implicit none 
     

    print*, '*************** MOLECULAR DYNAMICS **************'
    print*, 'STEP   TEMPERATURE    KINETIC ENERGY    POTENTIAL ENERGY'

    do step = 1, Nsteps
        call BOND_COND
        CALL MOVE_VERLET
       if (mod(step,100) .eq. 0 ) then   
           call temperature_control
           write(*,'(4X,I10,3F10.5)') step, temperature, en_kin, sum(en_pot(:))/npart
       endif
        if (mod(step,100) .eq. 0 .and. step .gt. Nsteps_equil) then 
		  CALL OVITO
		endif  
    enddo 

end subroutine dynamics