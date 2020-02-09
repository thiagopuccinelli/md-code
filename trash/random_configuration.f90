subroutine random_system 
    use variables_definition
    implicit none 
    integer j, iseed, i
    real ran2
    iseed = 918394
    do j = 1, npart 
        r(j,1) = ran2(iseed) * BOX
        r(j,2) = ran2(iseed) * BOX
        r(j,3) = ran2(iseed) * BOX
    enddo 
!    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! Step to generate the .xyz file for ovito
    open(200, file='random_config.dump')
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


end subroutine 