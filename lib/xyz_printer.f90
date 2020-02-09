subroutine OVITO
    use variables_definition  
    implicit none 
    integer :: i, k 

    open(10, file='sim.dump')
    write(10,'(A,3f6.2)') 'ITEM: TIMESTEP'
    write(10,'(I8)') step 
    write(10,'(A,3f6.2)') 'ITEM: NUMBER OF ATOMS'  
    write(10,'(I8)') npart
    write(10,'(A,3f6.2)') 'ITEM: BOX BOUNDS pp pp pp'
    write(10,'(2F25.10)') 0.d0, BOX
    write(10,'(2F25.10)') 0.d0, BOX
    write(10,'(2F25.10)') 0.d0, BOX
  	write(10,'(A,3f6.2)') 'ITEM: ATOMS id type x y z vx vy vz'
  	do i = 1, npart 
    	write(10,2000) i, '1' , r(i,:), v(i,:)
  	end do
2000 format(8X,I10,A3,6F25.10)
    return 
end subroutine