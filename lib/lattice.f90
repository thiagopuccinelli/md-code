subroutine lattice
    use variables_definition  
    implicit none 
    integer i, j, k, itel, n 
    real del, dx, dy, dz 
	n = int ( npart ** (1./3.) ) + 1 
	if (n .eq. 0) n = 1 
	del = BOX / DBLE (n)
	itel = 0 
	dx = - del 
	do i = 1, n 
		dx = dx + del 
		dy = - del 
		do j = 1, n 
			dy = dy + del 
			dz = - del 
			do k = 1, n 
				dz = dz + del 
				if (itel .lt. npart ) then 
					itel = itel + 1 
					r(itel,1) = dx + 0.5
					r(itel,2) = dy + 0.5
					r(itel,3) = dz + 0.5 
				endif 
			enddo
		enddo
	enddo
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! Step to generate the .xyz file for ovito
    open(200, file='lattice.dump')
    write(200,'(A,3f6.2)') 'ITEM: TIMESTEP'
    write(200,'(I8)') 0 
    write(200,'(A,3f6.2)') 'ITEM: NUMBER OF ATOMS'  
    write(200,'(I8)') npart
    write(200,'(A,3f6.2)') 'ITEM: BOX BOUNDS pp pp pp'
    write(200,'(2F25.10)') 0.d0, BOX
    write(200,'(2F25.10)') 0.d0, BOX
    write(200,'(2F25.10)') 0.d0, BOX
  	write(200,'(A,3f6.2)') 'ITEM: ATOMS id type x y z'
  	do i = 1, npart 
    	write(200,1000) i, '1' , r(i,:)
  	end do
1000 format(5X,I4,A3,3F25.10)
    close(200)

end subroutine lattice 