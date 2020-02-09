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