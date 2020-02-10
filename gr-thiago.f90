program gr 
   implicit none 

   integer, parameter :: N = 1000
   integer, parameter :: snaps = 100 
   real(8), parameter :: box_length = 11.0064241630
   integer, parameter :: Nhis = 2.**8. 
   real(8), parameter :: time = 1.0

   real(8) :: x(snaps,N), y(snaps,N), z(snaps,N), vx(snaps,N), vy(snaps,N), vz(snaps,N)
   integer :: type(snaps,N), id(snaps,N)
   integer i, j, k 
   integer dumb1,dumb2,dumb3
   real(8) :: avg(Nhis),r(Nhis)
   real(8) :: dr(snaps, N), x0,y0, z0, msd, msdvec(snaps), timevec(snaps)


   do j = 1, snaps 
      read(*,*) dumb1 
      read(*,*) dumb2 
      read(*,*) dumb3
      do i = 1, N 
         read(*,*) id(j,i), type(j,i), x(j,i), y(j,i), z(j,i), vx(j,i), vy(j,i), vz(j,i)
      enddo 
   enddo 

   dr = 0.0d0 
   do i = 2, snaps
      do j = 1, N 
         x0 = x(i-1,j)
         y0 = y(i-1,j)
         z0 = z(i-1,j)
         dr(i,j) = (x(i,j)-x0)**2.0d0+(y(i,j)-y0)**2.0d0+(z(i,j)-z0)**2.0d0
      enddo
   enddo 

   OPEN(unit=10,file='teste_read_MSD.dat',action="write")

   write(10,*) 0.0000001, 0.000001
   msdvec = 0.0d0
   timevec = 0.0d0
   msd = 0.0d0
   do j=2,snaps
      do i=1, N
         msd = msd + dr(j,i)
      enddo
      write(10,*) (j-1)*time , msd/N/2.0d0
      msdvec(j) = msd/N/2.0d0
      timevec(j) = (j-1)*time
   enddo


   call RADIAL_XYZ(N,snaps,x,y,z,r,avg,box_length,Nhis)

end program 

SUBROUTINE RADIAL_XYZ(Npart,Nframes,x,y,z,r,avg, LL, Nhis)
IMPLICIT NONE

INTEGER::Nhis,i,j,k,Nframes,ig,N1,sumN1,frames, Npart

DOUBLE PRECISION::rr,delg,lbox,pi,xr,yr,zr,r2,vb,nid,rho,LL,ddd
DOUBLE PRECISION,DIMENSION(Nframes,Npart)::x,y,z
INTEGER, DIMENSION(Nframes,Npart):: type
DOUBLE PRECISION,DIMENSION(10000, Nhis)::gr
DOUBLE PRECISION,DIMENSION(Nhis),INTENT(OUT)::r,avg

delg=LL/(2.*Nhis)
pi=4*ATAN(1.)
N1 = Npart
rho = N1/(LL**3.0d0)
gr = 0.0d0
avg(:)=0.d0

DO k=1,Nframes
   DO i=1,Npart-1
      DO j=i+1,Npart
         xr=x(k,i)-x(k,j)
         yr=y(k,i)-y(k,j)
         zr=z(k,i)-z(k,j)
            
         xr=xr-LL*(NINT(xr/LL))
         yr=yr-LL*(NINT(yr/LL))
         zr=zr-LL*(NINT(zr/LL))
         r2=xr*xr+yr*yr+zr*zr
         rr=SQRT(r2)

         IF(rr.LT.LL/2.d0)THEN
            ig=ceiling(rr/delg)
            gr(k,ig)=gr(k,ig)+2.
         END IF
      END DO
   END DO
END DO




DO j=1,Nhis
   DO i=1,Nframes
      r(j)=delg*(j+0.5)
      vb=((j+1)**3.-j**3.)*delg**3.
      nid=(4./3.)*pi*vb*rho
      gr(i,j)=gr(i,j)/(N1*nid)
   END DO
END DO





DO i=1,Nhis
   DO j=1,Nframes
      avg(i)=avg(i)+gr(j,i)
   END DO
END DO

OPEN(unit=2,file='teste_read_xyz.dat',action="write")

DO i=1,Nhis
   WRITE(2,'(2(f17.10,1X))')r(i),avg(i)/Nframes
END DO


END SUBROUTINE RADIAL_XYZ
