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
    real, parameter :: gamma = 10.0d0 
    integer, parameter :: dim = 3  
    real :: r(npart,dim), f(npart,dim), v(npart,dim)
    real :: en_pot(npart)
end module vars
program md 
    use vars 
    implicit none 
    integer :: step
    real :: temp, en_kin 
    r = 0 
    f = 0 
    v = 0 
    en_pot = 0.d0 
    en_kin = 0.d0 
    call INIT
    call FORCE
    print*, '*************** MOLECULAR DYNAMICS **************'
    print*, 'STEP   TEMPERATURE    KINETIC ENERGY    POTENTIAL ENERGY'
    do step = 1, Nsteps
        call BOND_COND
        CALL MOVE_VERLET(step)
       if (mod(step,100) .eq. 0 ) then   
           call temperature(step,temp, en_kin)
           write(*,'(4X,I10,3F10.5)') step, temp, en_kin, sum(en_pot(:))/npart
       endif
        if (mod(step,100) .eq. 0 .and. step .gt. Nsteps_equil) then 
		    CALL OVITO(step)
		endif  
    enddo 
end program 

subroutine INIT
    use vars 
    implicit none 
    integer :: i, k 
    integer :: iseed
    real vel
    real fs, v2 !for the setting zero center of mass
    real, parameter :: temp = 1.d0 
    real gasdev 
    iseed = 387489
    call LATTICE
    do i = 1, npart 
        do k = 1, dim 
            vel = gasdev(iseed)
            v(i,k) = vel
        enddo 
    enddo 
    return 
end subroutine INIT 
subroutine LATTICE
	use vars
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
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
end subroutine LATTICE

subroutine FORCE
    use vars 
    implicit none 
    integer :: i, j,k 
    real :: virial 
    real dev1, dev, g, gasdev, phi 
    real rij(dim),sij(dim)
    real r2i, r6i, ff, r2, r12i 
    real, dimension(dim) :: LANGEVIN
    integer ::  iseed 
    iseed = 9381283
    !CALL RANDOM_SEED()
    f = 0.d0
    en_pot = 0.d0 

    DO i=1,npart-1
        DO j=i+1,npart
        sij(1:dim) = r(i, 1:dim) - r(j, 1:dim)
        sij(1:dim) = sij(1:dim) - BOX * (NINT(sij(1:dim) / BOX))
        r2 = dot_product (sij(1:dim), sij(1:dim))
         IF(r2.LT.Rcutoff**2)THEN
            r2i = 1 / r2
            r6i = r2i ** 3. 
            r12i = r6i ** 2. 
            phi = 4.d0 * (r12i - r6i) - phicutoff
            ff = 24.d0 * r2i * (2.d0 * r12i - r6i)
            en_pot(i) = en_pot(i) + 0.5d0 * phi  
            en_pot(j) = en_pot(j) + 0.5d0 * phi  
            f(i,1:dim) = f(i,1:dim) + ff * sij(1:dim)
            f(j,1:dim) = f(j,1:dim) - ff * sij(1:dim)
         END IF
      END DO
   END DO


    do i = 1, npart 
        dev1 = 6. * temprqs / dt 
        dev = sqrt(dev1 * gamma)
        g = gasdev(iseed)
        g = g * dev 
        call genUnitRandVector(LANGEVIN,iseed)
        f(i,1:dim) = f(i,1:dim) + LANGEVIN(1:dim) * g
		f(i,1:dim) = f(i,1:dim) - gamma * v(i,1:dim)
	enddo 
    return 
end subroutine FORCE 

subroutine BOND_COND
    use vars 
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

subroutine MOVE_VERLET(step)
! 	** VELOCITY VERLET ALGORITHM **
!	** AT THE START OF A TIMESTEP, MOVE_VELVERLET_A IS CALLED TO ADVANCE THE 
!	** THE POSITIONS AND 'HALF-ADVANCE' THE VELOCITIES. THEN THE FORCE ROUTINE 
!	** IS CALLED, AND THIS IS FOLLOWED BY MOVE_VELVERLET_B WHICH COMPLETES THE 
!	** ADVANCEMENT OF VELOCITIES 
    use vars 
    implicit none 
    integer :: step, i, k
    real :: dt2, dt22, time
    
    dt2 = dt / 2.0 
    dt22 = dt * dt2 
    time = step * dt 

    do k = 1, dim 
        do i = 1, npart
            r(i,k) = r(i,k) + v(i,k) * dt + dt22 * (f(i,k)) / mass 
        enddo  
        call FORCE
        do i = 1, npart 
            v(i,k) =  v(i,k) + dt2 * f(i,k) / mass
        enddo 
    enddo 
    return 
end subroutine MOVE_VERLET

subroutine OVITO(step)
    use vars 
    implicit none 
    integer :: step, i, k 

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

subroutine temperature(step,temp, en_kin)
    use vars 
    implicit none 
    integer :: i, step, k  
    real :: temp, v2, en_kin
    v2 = 0.d0 
    temp = 0.d0
    en_kin = 0.d0
    do k = 1, dim  
        do i = 1, npart
            v2 = v2 + v(i,k) ** 2
        enddo 
    enddo 
    temp  = v2 / (3*npart)
    v2 = v2 / npart 
    en_kin = 0.5d0 * mass * v2 
    open(30, file='ENERGY.dat')
    write(30, 3000) step, temp, en_kin, sum(en_pot(:))/npart
3000 format(4X, I6, 3F20.5)
    return 
end subroutine temperature


!--------------------------------------------------------------------------------
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
  