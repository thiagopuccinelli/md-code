module vars 
    implicit none 
    integer, parameter :: Nsteps = 15000
    integer, parameter :: Nsteps_equil = 5000
    real, parameter :: dt = 0.001d0 
    real, parameter :: t = dt * Nsteps 
    integer, parameter :: npart = 1000 
    real, parameter :: mass = 1.0d0
    real, parameter :: eps = 1.0d0
    real, parameter :: sigma = 1.0d0
    real, parameter :: rho = 0.85d0
    real, parameter :: BOX = (npart/rho) ** (1./3.)
    real, parameter :: Rcutoff = 2.5d0   ! cutoff distance 2.5*sigma
  	real, parameter :: phicutoff =  1./(Rcutoff**12) - 1.d0/(Rcutoff**6)    ! potential at cutoff
    real, parameter :: temprqs = 1.d0
    real, parameter :: gamma = 1.0d0  
    integer xx(npart), yy(npart), zz(npart) !for unfolding positions 
!variables for MSD measure:
	integer, parameter :: TDIFMAX = Nsteps - Nsteps_equil 
	integer, parameter :: timemax = Nsteps - Nsteps_equil  
	integer, parameter :: it0 =	10
    integer, parameter :: t0max = (Nsteps - Nsteps_equil)/10
! variables for the neighbouring listing 
	integer, parameter :: M = floor(BOX/Rcutoff)
	integer, parameter :: NCELL = M * M * M 
	integer, parameter :: MAPSIZ = 13 * NCELL 
	integer LIST(50000), HEAD(NCELL), MAP(MAPSIZ)
end module vars
program md 
    use vars 
    implicit none 
    integer :: step
    real, dimension(npart) :: RX, RY, RZ
    real, dimension(npart) :: VX, VY, VZ
    real, dimension(npart) :: FX, FY, FZ
    real :: temp, en_pot, en_kin 
    !	** VARIABLES FOR THE DIFFUSION CALCULATION ** 
	integer, parameter :: nsamp = 1
    RX = 0.d0 
    RY = 0.d0
    RZ = 0.d0
    VX = 0.d0
    VY = 0.d0
    VZ = 0.d0
    FX = 0.d0
    FY = 0.d0
    FZ = 0.d0
    en_pot = 0.d0 
    en_kin = 0.d0 
    call INIT(RX,RY,RZ,VX,VY,VZ) 
    call MAPS 
    call FORCE(RX,RY,RZ,VX,VY,VZ,FX,FY,FZ,en_pot)
    call dif(0, nsamp, RX, RY, RZ)
    print*, '*************** MOLECULAR DYNAMICS **************'
    print*, 'STEP   TEMPERATURE    KINETIC ENERGY    POTENTIAL ENERGY'
    do step = 1, Nsteps
        call BOND_COND(RX,RY,RZ)
        call LINKS(RX,RY,RZ) 
        CALL MOVE_VERLET(RX,RY,RZ,VX,VY,VZ,FX,FY,FZ,step,en_pot)
        if (mod(step,100) .eq. 0 ) then   
            call temperature(step,temp, VX, VY, VZ, en_pot, en_kin)
            write(*,'(4X,I10,3F10.5)') step, temp, en_kin, en_pot 
        endif
        if (mod(step,100) .eq. 0 .and. step .gt. Nsteps_equil) then 
			CALL OVITO(step,RX,RY,RZ,VX,VY,VZ)
		endif  
        !	** DIFFUSION CALCULATION ** 
		if (mod(step,nsamp) .eq. 0 .and. step .gt. Nsteps_equil) then 
			call dif(1, nsamp, RX, RY, RZ)
		endif
    enddo 
    !	** PRINT DIFFUSION ANALYSIS ** 
	call dif(2, nsamp, RX, RY, RZ) 
end program 

subroutine INIT(RX,RY,RZ,VX,VY,VZ)
    use vars 
    implicit none 
    integer :: i
    integer :: iseed
    real sumvx, sumvy, sumvz, sumvx2, sumvy2, sumvz2
    real velx, vely, velz
    real fs, v2 !for the setting zero center of mass
    real, dimension(npart) :: RX, RY, RZ
    real, dimension(npart) :: VX, VY, VZ
    real, parameter :: temp = 1.d0 
    real gasdev 
    iseed = 387489
    sumvx = 0.d0
	sumvy = 0.d0 
	sumvz = 0.d0 
	sumvx2 = 0.d0 
	sumvy2 = 0.d0 
	sumvz2 = 0.d0 
    RX = 0.d0 
    RY = 0.d0
    RZ = 0.d0
    VX = 0.d0
    VY = 0.d0
    VZ = 0.d0
    call LATTICE(RX,RY,RZ)
    do i = 1, npart 
        velx = gasdev(iseed)
        vely = gasdev(iseed)
        velz = gasdev(iseed)
        VX(i) = velx 
        VY(i) = vely 
        VZ(i) = velz 
    ! compute the center of mass velocity 
        sumvx = sumvx + VX(i)
        sumvy = sumvy + VY(i)
        sumvz = sumvz + VZ(i)
        sumvx2 = sumvx2 + VX(i) ** 2.d0 
        sumvy2 = sumvy2 + VY(i) ** 2.d0 
        sumvz2 = sumvz2 + VZ(i) ** 2.d0 
    enddo 
! CENTER OF MASS 
    sumvx = sumvx / npart 
    sumvy = sumvy / npart 
    sumvz = sumvz / npart 
! MEAN QUADRATIC VELOCITY 
    sumvx2 = sumvx2 / npart
    sumvy2 = sumvy2 / npart
    sumvz2 = sumvz2 / npart
    v2 = sumvx2 + sumvy2 + sumvz2 
    fs = sqrt( 3 * temp / v2 ) ! Scaling velocities factor
!	** Fixing the center of mass **
    do i = 1, npart 
        VX(i) = ( VX(i) - sumvx ) * fs 
        VY(i) = ( VY(i) - sumvy ) * fs 
        VZ(i) = ( VZ(i) - sumvz ) * fs
        RX(i) = RX(i) - VX(i) * dt
        RY(i) = RY(i) - VY(i) * dt
        RZ(i) = RZ(i) - VZ(i) * dt  
    enddo 
    return 
end subroutine INIT 
subroutine LATTICE(RX,RY,RZ)
	use vars
	implicit none 
	integer i, j, k, itel, n 
    real del, dx, dy, dz 
    real, dimension(npart) :: RX, RY, RZ
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
					RX(itel) = dx + 0.5
					RY(itel) = dy + 0.5
					RZ(itel) = dz + 0.5 
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
    	write(200,1000) i, '1' , RX(i), RY(i), RZ(i)
  	end do
1000 format(5X,I4,A3,3F25.10)
    close(200)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
end subroutine LATTICE

subroutine FORCE(RX,RY,RZ,VX,VY,VZ,FX,FY,FZ,en_pot)
    use vars 
    implicit none 
    integer :: i, j, ICELL, JCELL0, JCELL, NABOR
    real, dimension(npart) :: RX, RY, RZ
    real, dimension(npart) :: VX, VY, VZ
    real, dimension(npart) :: FX, FY, FZ
    real RXI, RYI, RZI, FXI, FYI, FZI
    real RXIJ, RYIJ, RZIJ, RIJSQ, SR2, SR6, SR12
    real V, VIJ, W, WIJ, FXIJ, FYIJ, FZIJ, FIJ
    real :: en_pot, virial 
    real dev1, dev, g, gasdev
    real, dimension(3) :: F  
    integer ::  iseed 
    iseed = 9381283
    !CALL RANDOM_SEED()
    do i = 1, npart 
        FX(i) = 0.d0
        FY(i) = 0.d0
        FZ(i) = 0.d0
    enddo 
    V = 0.d0 
    W = 0.d0 
    en_pot = 0.d0 

!	** LOOP OVER ALL CELLS BEGINS ** 
    do 5000 ICELL = 1, NCELL
        I = HEAD(ICELL)

!	** LOOP OVER ALL MOLECULES IN THE CELL ** 
1000		if (I .GT. 0) then  
				RXI = RX(i)
				RYI = RY(i)
				RZI = RZ(i)
				FXI = FX(i)
				FYI = FY(i)
				FZI = FZ(i)
!		** LOOP OVER ALL MOLECULES BELOW I IN THE CURRENT CELL **
				J = LIST(I)
2000			IF (J .GT. 0) then  
					RXIJ = RXI - RX(J)
					RYIJ = RYI - RY(J)
					RZIJ = RZI - RZ(J)
					RXIJ = RXIJ - BOX * nint( RXIJ / BOX )
					RYIJ = RYIJ - BOX * nint( RYIJ / BOX )
					RZIJ = RZIJ - BOX * nint( RZIJ / BOX )
					RIJSQ = RXIJ ** 2. + RYIJ ** 2. + RZIJ ** 2. 

					if (RIJSQ .lt. Rcutoff**2) then 
						SR2 = 1 / RIJSQ 
						SR6 = SR2 * SR2 * SR2
						SR12 = SR6 * SR6
						VIJ = SR12 - SR6 
						V = V + VIJ !- phicutoff
						WIJ = VIJ + SR12
						W = W + WIJ 
						FIJ = WIJ / RIJSQ
						FXIJ = FIJ * RXIJ
						FYIJ = FIJ * RYIJ
						FZIJ = FIJ * RZIJ 
						FXI = FXI + FXIJ 
						FYI = FYI + FYIJ
						FZI = FZI + FZIJ
						FX(J) = FX(J) - FXIJ 
						FY(J) = FY(J) - FYIJ 
						FZ(J) = FZ(J) - FZIJ 
					endif 
					J = LIST(J)
					go to 2000 
			endif 
!		** LOOP OVER NEIGHBOURING CELLS ** 
			JCELL0 = 13 * (ICELL - 1)
			do 4000 NABOR = 1, 13 
				JCELL = MAP ( JCELL0 + NABOR )
!		** LOOP OVER ALL MOLECULES IN NEIGHBOURING CELLS **
				J = HEAD(JCELL)

3000			if (J .NE. 0) then 
					RXIJ = RXI - RX(J)
					RYIJ = RYI - RY(J)
					RZIJ = RZI - RZ(J)
					RXIJ = RXIJ - BOX * nint( RXIJ / BOX )
					RYIJ = RYIJ - BOX * nint( RYIJ / BOX )
					RZIJ = RZIJ - BOX * nint( RZIJ / BOX )
					RIJSQ = RXIJ ** 2. + RYIJ ** 2. + RZIJ ** 2.
						if (RIJSQ .lt. Rcutoff**2) then 
							SR2 = 1 / RIJSQ 
							SR6 = SR2 * SR2 * SR2
							SR12 = SR6 * SR6
							VIJ = SR12 - SR6 
							V = V + VIJ !- phicutoff
							WIJ = VIJ + SR12
							W = W + WIJ 
							FIJ = WIJ / RIJSQ
							FXIJ = FIJ * RXIJ
							FYIJ = FIJ * RYIJ
							FZIJ = FIJ * RZIJ 
							FXI = FXI + FXIJ 
							FYI = FYI + FYIJ
							FZI = FZI + FZIJ
							FX(J) = FX(J) - FXIJ 
							FY(J) = FY(J) - FYIJ 
							FZ(J) = FZ(J) - FZIJ 
						endif 
						J = LIST(J)
						go to 3000
				endif 
4000		continue
			
!	** INNER LOOP ENDS **
	FX(i) = FXI 
	FY(i) = FYI 
	FZ(i) = FZI 

	I = LIST(I) 

	go to 1000
	endif 
5000 continue
    do i = 1, npart 
        FX(i) = FX(i) * eps * 24
        FY(i) = FY(i) * eps * 24
        FZ(i) = FZ(i) * eps * 24
    enddo 
    do i = 1, npart 
        dev1 = 6. * temprqs / dt 
        dev = sqrt(dev1 * gamma)
        g = gasdev(iseed)
        g = g * dev 
        call genUnitRandVector(F,iseed)
        FX(I) = FX(I) + F(1) * g
		FY(I) = FY(I) + F(2) * g 
		FZ(I) = FZ(I) + F(3) * g 
		FX(I) = FX(I) - gamma * VX(I)
		FY(I) = FY(I) - gamma * VY(I)
		FZ(I) = FZ(I) - gamma * VZ(I)
	enddo 
    en_pot = V * 4. * eps - 4. * eps * phicutoff
    en_pot = en_pot / npart 
    W = W * eps * 24. / 3.0 
    virial = -W
    return 
end subroutine FORCE 

subroutine BOND_COND(RX,RY,RZ)
    use vars 
    real, dimension(npart) :: RX,RY,RZ
    integer :: i 

    do i = 1, npart
        if (RX(i) .lt. 0.d0) then 
            RX(i) = RX(i) + BOX 
            xx(i) = xx(i) - 1 
        endif  
        if (RY(i) .lt. 0.d0) then 
              RY(i) = RY(i) + BOX
              yy(i) = yy(i) - 1 
        endif 
        if (RZ(i) .lt. 0.d0) then 
              RZ(i) = RZ(i) + BOX
              zz(i) = zz(i) - 1 
        endif 
    end do
    
    do i = 1, npart
         if (RX(i) .gt. BOX) then 
             RX(i) = RX(i) - BOX
             xx(i) = xx(i) + 1 
         endif 
         if (RY(i) .gt. BOX) then 
             RY(i) = RY(i) - BOX
             yy(i) = yy(i) + 1 
         endif 
         if (RZ(i) .gt. BOX) then 
             RZ(i) = RZ(i) - BOX
             zz(i) = zz(i) + 1 
         endif 
    end do  
    return 
end subroutine BOND_COND

subroutine MOVE_VERLET(RX,RY,RZ,VX,VY,VZ,FX,FY,FZ,step,en_pot)
! 	** VELOCITY VERLET ALGORITHM **
!	** AT THE START OF A TIMESTEP, MOVE_VELVERLET_A IS CALLED TO ADVANCE THE 
!	** THE POSITIONS AND 'HALF-ADVANCE' THE VELOCITIES. THEN THE FORCE ROUTINE 
!	** IS CALLED, AND THIS IS FOLLOWED BY MOVE_VELVERLET_B WHICH COMPLETES THE 
!	** ADVANCEMENT OF VELOCITIES 
    use vars 
    implicit none 
    integer :: step, i
    real, dimension(npart) :: RX,RY,RZ
    real, dimension(npart) :: VX,VY,VZ
    real, dimension(npart) :: FX,FY,FZ
    real :: dt2, dt22, time, en_pot 
    
    dt2 = dt / 2.0 
    dt22 = dt * dt2 
    time = step * dt 

    do i = 1, npart
        RX(i) = RX(i) + VX(i) * dt + dt22 * (FX(i)) / mass 
        RY(i) = RY(i) + VY(i) * dt + dt22 * (FY(i)) / mass 
        RZ(i) = RZ(i) + VZ(i) * dt + dt22 * (FZ(i)) / mass 
        VX(i) = VX(i) + dt2 * FX(i) / mass 
        VY(i) = VY(i) + dt2 * FY(i) / mass 
        VZ(i) = VZ(i) + dt2 * FZ(i) / mass 
    enddo  
    call FORCE(RX,RY,RZ,VX,VY,VZ,FX,FY,FZ,en_pot)
    do i = 1, npart 
        VX(i) =  VX(i) + dt2 * FX(i) / mass
        VY(i) =  VY(i) + dt2 * FY(i) / mass
        VZ(i) =  VZ(i) + dt2 * FZ(i) / mass
    enddo 
    return 
end subroutine MOVE_VERLET

subroutine OVITO(step,RX,RY,RZ,VX,VY,VZ)
    use vars 
    implicit none 
    integer :: step, i 
    real, dimension(npart) :: RX,RY,RZ
    real, dimension(npart) :: VX,VY,VZ

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
    	write(10,2000) i, '1' , RX(i), RY(i), RZ(i), VX(i), VY(i), VZ(i)
  	end do
2000 format(8X,I10,A3,6F25.10)
    return 
end subroutine

subroutine temperature(step,temp, VX, VY, VZ, en_pot, en_kin)
    use vars 
    implicit none 
    integer :: i, step 
    real :: temp, v2, en_kin, en_pot
    real, dimension(npart) :: VX, VY, VZ
    v2 = 0.d0 
    temp = 0.d0
    en_kin = 0.d0 
    do i = 1, npart
		v2 = v2 + VX(i) ** 2. + VY(i) ** 2. + VZ(i) ** 2. 
    enddo 
    temp  = v2 / (3*npart)
    v2 = v2 / npart 
    en_kin = 0.5d0 * mass * v2 
    open(30, file='ENERGY.dat')
    write(30, 3000) step, temp, en_kin, en_pot
3000 format(4X, I6, 3F20.5)
    return 
end subroutine temperature

!********************************************************************************
!** FICHE F.20.  ROUTINES TO CONSTRUCT AND USE CELL LINKED-LIST METHOD         **
!** This FORTRAN code is intended to illustrate points made in the text.       **
!** To our knowledge it works correctly.  However it is the responsibility of  **
!** the user to test it, if it is to be used in a research application.        **
!********************************************************************************

!    *******************************************************************
!    ** CONSTRUCTION OF CELL LINKED-LISTS AND USE IN FORCE ROUTINE.   **
!    **                                                               **
!    ** REFERENCES:                                                   **
!    **                                                               **
!    ** QUENTREC AND BROT, J. COMPUT. PHYS. 13, 430, 1975.            **
!    ** HOCKNEY AND EASTWOOD, COMPUTER SIMULATION USING PARTICLES,    **
!    **    MCGRAW HILL, 1981.                                         **
!    **                                                               **
!    ** ROUTINES SUPPLIED:                                            **
!    **                                                               **
!    ** SUBROUTINE MAPS                                               **
!    **    SETS UP MAP OF CELL STRUCTURE FOR USE IN FORCE             **
!    ** SUBROUTINE LINKS ( RCUT )                                     **
!    **    SETS UP HEAD OF CHAIN ARRAY AND LINKED LIST                **
!    ** SUBROUTINE FORCE ( SIGMA, RCUT, V, W )                        **
!    **    CALCULATES FORCES USING A LINKED LIST                      **
!    **                                                               **
!    ** USAGE:                                                        **
!    **                                                               **
!    ** SUBROUTINE MAPS IS CALLED ONCE AT THE START OF A SIMULATION   **
!    ** TO ESTABLISH CELL NEIGHBOUR IDENTITIES.  AT EACH TIMESTEP,    **
!    ** SUBROUTINE LINKS IS CALLED TO SET UP THE LINKED LIST AND THIS **
!    ** IS IMMEDIATELY USED BY SUBROUTINE FORCE.                      **
!    *******************************************************************

subroutine MAPS 
	use vars 
	implicit none 

!   *******************************************************************
!    ** ROUTINE TO SET UP A LIST OF NEIGHBOURING CELLS                **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER M                  NUMBER OF CELLS IN EACH DIRECTION  **
!    ** INTEGER MAPSIZ             SIZE OF CELL-CELL MAP              **
!    ** INTEGER MAP(MAPSIZ)        LIST OF NEIGHBOURING CELLS         **
!    **                                                               **
!    ** USAGE:                                                        **
!    **                                                               **
!    ** THIS SUBROUTINE SETS UP A LIST OF THE THIRTEEN NEIGHBOURING   **
!    ** CELLS OF EACH OF THE SMALL CELLS IN THE CENTRAL BOX. THE      **
!    ** EFFECTS OF THE PERIODIC BOUNDARY CONDITIONS ARE INCLUDED.     **
!    ** THE SUBROUTINE IS CALLED ONCE AT THE BEGINNING OF THE         **
!    ** SIMULATION AND THE MAP IS USED IN THE FORCE SUBROUTINE        **
!    *******************************************************************

	integer IX, IY, IZ, IMAP, ICELL 

!	** STATEMENT TO GIVE CELL INDEX ** 
	ICELL(IX, IY, IZ) = 1 + MOD (IX - 1 + M, M) + MOD(IY - 1 + M, M) * M &
	                         + MOD (IZ - 1 + M, M) * M * M 

!	** FIND HALF THE NEAREST NEIGHBOURS OF EACH CELL **
	do 50 IZ = 1, M 
		do 40 IY = 1, M 
			do 30 IX = 1, M 
				IMAP = ( ICELL (IX, IY, IZ) - 1 ) * 13 

				MAP( IMAP + 1  ) = ICELL( IX + 1, IY    , IZ     )
                MAP( IMAP + 2  ) = ICELL( IX + 1, IY + 1, IZ     )
                MAP( IMAP + 3  ) = ICELL( IX    , IY + 1, IZ     )
                MAP( IMAP + 4  ) = ICELL( IX - 1, IY + 1, IZ     )
                MAP( IMAP + 5  ) = ICELL( IX + 1, IY    , IZ - 1 )
                MAP( IMAP + 6  ) = ICELL( IX + 1, IY + 1, IZ - 1 )
                MAP( IMAP + 7  ) = ICELL( IX    , IY + 1, IZ - 1 )
                MAP( IMAP + 8  ) = ICELL( IX - 1, IY + 1, IZ - 1 )
                MAP( IMAP + 9  ) = ICELL( IX + 1, IY    , IZ + 1 )
                MAP( IMAP + 10 ) = ICELL( IX + 1, IY + 1, IZ + 1 )
                MAP( IMAP + 11 ) = ICELL( IX    , IY + 1, IZ + 1 )
                MAP( IMAP + 12 ) = ICELL( IX - 1, IY + 1, IZ + 1 )
                MAP( IMAP + 13 ) = ICELL( IX    , IY    , IZ + 1 )
30          continue 

40      continue 

50  continue 

	RETURN 

end subroutine MAPS

subroutine LINKS(RX,RY,RZ)
	use vars
    implicit none 
    real, dimension(npart) :: RX,RY,RZ



!    *******************************************************************
!    ** ROUTINE TO SET UP LINKED LIST AND THE HEAD OF CHAIN ARRAYS    **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER N                  NUMBER OF ATOMS                    **
!    ** INTEGER M                  NUMBER OF CELLS IN EACH DIRECTION  **
!    ** INTEGER NCELL              TOTAL NUMBER OF CELLS (M**3)       **
!    ** INTEGER LIST(N)            LINKED LIST OF ATOMS               **
!    ** INTEGER HEAD(NCELL)        HEAD OF CHAIN FOR EACH CELL        **
!    ** REAL    RX(N),RY(N),RZ(N)  POSITIONS                          **
!    ** REAL    RCUT               THE CUTOFF DISTANCE FOR THE FORCE  **
!    **                                                               **
!    ** USAGE:                                                        **
!    **                                                               **
!    ** EACH ATOM IS SORTED INTO ONE OF THE M**3 SMALL CELLS.         **
!    ** THE FIRST ATOM IN EACH CELL IS PLACED IN THE HEAD ARRAY.      **
!    ** SUBSEQUENT ATOMS ARE PLACED IN THE LINKED LIST ARRAY.         **
!    ** ATOM COORDINATES ARE ASSUMED TO BE BETWEEN -0.5 AND +0.5.     **
!    ** THE ROUTINE IS CALLED EVERY TIMESTEP BEFORE THE FORCE ROUTINE.**
!    *******************************************************************


	real CELLI, CELL
	integer ICELL, I 


!	** ZERO HEAD OF CHAIN ARRAY 
	do ICELL = 1, NCELL 
		HEAD(ICELL) = 0 
	enddo 

	CELLI = DBLE (M)
	CELL = BOX / CELLI 

	if (CELL .LT. Rcutoff) then 
		stop 'CELL SIZE TOO SMALL FOR CUTOFF'
	endif 

!	** SORT ALL ATOMS ** 
	do I = 1, npart 
		ICELL = 1 + INT( ( RX(I) / BOX ) * CELLI ) & 
				  + INT( ( RY(I) / BOX ) * CELLI ) * M &
				  + INT( ( RZ(I) / BOX ) * CELLI ) * M * M 

		LIST(I) = HEAD(ICELL)
		HEAD(ICELL) = I 
	enddo 

	RETURN 
end subroutine LINKS 
!--------------------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine dif(switch, nsamp, RX, RY, RZ)
	use vars
    implicit none 
    integer :: switch, nsamp, delt, i, tstep 
    real  r2a, r2asum 
    integer tvacf, t0, tt0, ngr
	real :: dtime 
	real time 
	integer, dimension(t0max) :: time0 
	integer ntime(timemax)
	real X0(npart,t0max), Y0(npart,t0max), Z0(npart,t0max), r2t(timemax)
	real :: RXU(npart), RYU(npart), RZU(npart)
    real, dimension(npart) :: RX, RY, RZ
    save tvacf, t0, tt0, ngr 
	if (switch .eq. 0) then 
        tvacf = 0 
		delt = 0 
		dtime = dt * DBLE(nsamp) 
		do i = 1, timemax 
			ntime(i) = 0
			r2t(i) = 0 
		enddo 
	else if (switch .eq. 1) then 
		tvacf = tvacf + 1 
		dtime = dt * nsamp 
		if (mod(tvacf, it0) .eq. 0) then 
			t0 = t0 + 1 !update t = 0 
			tt0 = mod(t0 - 1, t0max) + 1 
			time0(tt0) = tvacf       !store the time t=0
			do i = 1, npart 
				X0(i, tt0) = RX(i) + xx(i) * BOX 
				Y0(i, tt0) = RY(i) + yy(i) * BOX 
				Z0(i, tt0) = RZ(i) + zz(i) * BOX 
			enddo
		endif 
		do tstep = 1, min(t0, t0max)
			delt = tvacf - time0(tstep) + 1 
			if (delt.lt.timemax .and. dt*dtime .le. TDIFMAX) then 
				ntime(delt) = ntime(delt) + 1 
				r2a = 0.d0 
				r2asum = 0.d0 
				do i = 1, npart 
!				** UNFOLD POSITIONS ** 
					RXU(i) = RX(i) + xx(i) * BOX 
					RYU(i) = RY(i) + yy(i) * BOX 
					RZU(i) = RZ(i) + zz(i) * BOX 
					r2a = ( RXU(i) - X0(i,tstep)) ** 2. + ( RYU(i) - Y0(i,tstep)) ** 2. + ( RZU(i) - Z0(i,tstep)) ** 2.
					r2t(delt) = r2t(delt) + r2a 
					r2asum = r2asum + r2a 
				enddo
			endif 
		enddo 
	else if (switch .eq. 2) then 
		open(23, file='msd.dat') 
  		dtime = dt * nsamp 
  		do i = 1, timemax 
  			time = dtime * (i + 0.5)
  			if (ntime(i) .ne. 0) then 
  				r2t(i) = r2t(i) / ( npart * ntime(i) )
  				write(23, '(2F20.10)' ) time, r2t(i)
  			endif
  		enddo 
	endif 
	return 
end subroutine dif 
!---------------------------------------------------------------!
real function gasdev(idum)
  ! (C) Copr. 1986-92 Numerical Recipes Software .
  !  USES ran2
  integer idum
  integer iset
  real fac,gset,rsq,v1,v2,ran2
  SAVE iset,gset
  DATA iset/0/
  if (iset.eq.0) then
1    v1 = 2.*ran2(idum) - 1.
     v2 = 2.*ran2(idum) - 1.
     rsq = v1**2 + v2**2
     if(rsq.ge.1..or.rsq.eq.0.)goto 1
     fac = sqrt(-2.*log(rsq)/rsq)
     gset = v1*fac
     gasdev = v2*fac
     iset = 1
  else
     gasdev = gset
     iset = 0
  endif
  return
end function gasdev
subroutine genUnitRandVector(p, pSeed)   !average success rate is pi/4

    implicit none
  
    real p(3)
    real x, y, s, ran2
    integer pSeed
  
    s = 2.0
    do while ( s > 1.0)
       x = 2.0 * ran2(pSeed) - 1.0
       y = 2.0 * ran2(pSeed) - 1.0
       s = x*x + y*y
    enddo
    p(3) = 1.0 - 2.0 * s
    s = 2.0 * sqrt(1.0 - s)
    p(1) = s * x
    p(2) = s * y
  
    return
  
  end subroutine genUnitRandVector

   FUNCTION ran2(idum)

    implicit none
  
    INTEGER idum,idum2
    INTEGER IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
    real ran2,AM,EPS,RNMX
    PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&
         IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
         NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
    INTEGER j,k,iv(NTAB),iy
    SAVE iv,iy,idum2
    DATA idum2/123456789/, iv/NTAB*0/, iy/0/
    if (idum.le.0) then
       idum=max(-idum,1)
       idum2=idum
       do j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
       enddo
       iy=iv(1)
    endif
    k=idum/IQ1
    idum=IA1*(idum-k*IQ1)-k*IR1
    if (idum.lt.0) idum=idum+IM1
    k=idum2/IQ2
    idum2=IA2*(idum2-k*IQ2)-k*IR2
    if (idum2.lt.0) idum2=idum2+IM2
    j=1+iy/NDIV
    iy=iv(j)-idum2
    iv(j)=idum
    if(iy.lt.1)iy=iy+IMM1
    ran2=min(AM*iy,RNMX)
    return
  END FUNCTION ran2
  