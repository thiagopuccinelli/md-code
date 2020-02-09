subroutine FORCE
    use variables_definition  
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
