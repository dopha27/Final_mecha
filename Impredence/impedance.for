C	ImpedanceWeight	
      USE MSFLIB
C	double complex Z(200),zn(200),z1,ZD1,ZD2,ZD3,Z0,zeta,theta,s
	double complex Z(200),zn(200),z1,ZD1,ZD2,ZD3,Z0,zeta,theta,s
	complex table,dummy
	COMMON/TEMPO/table(41,41)
	dimension QR(100),QU(100)
	real*8 pi,aa,bb,eta
	Rroot=0.05
	dt=0.05
	N=20
	eps=1e-10
	eta=eps**(1.d0/2/N)
c	write(*,*) eps,eta
	pi=datan(1.d0)*4.d0
	alpa=0.9
	beta=0.6
	dummy=(-1,0.)
	pause
c	pi=datan(1.d0)*4.d0

c	m=5
        OPEN(5,FILE='cmplx1.dat',STATUS='UNKNOWN')

	do 20 m=0,2*N-1
	do 10 i=1,40
	do 10 j=1,10
		table(i,j)=dummy
c	write(5,*) dummy, table(i,j)
10	continue
	the=pi*m/N
	zeta=eta*cmplx(cos(the),sin(the),8)
	theta=zeta*zeta/2.-2*zeta+1.5
c	write(*,*) zeta,theta
	s=theta/dt
 	call impedance(s,Rroot,1,1,Z0,dummy)
c	Z(m)=Z0
	Z(m+1)=Z0
20    continue
	close(5)
c        OPEN(6,FILE='cmplx.dat',STATUS='UNKNOWN')
c	do 22 i=1,40
c	do 22 j=1,10
c	write(6,*) table(i,j)	  
c22	continue
 
 	do 30 i=0,N
		z1=(0.,0.)
	  do 25 j=0,2*N-1
		the=-i*j*pi/N
25	    z1=z1+Z(j+1)*cmplx(cos(the),sin(the),8)
	  zn(i+1)=z1*eta**(-i)/2/N
	write(*,*) i, zn(i+1)
30	continue
      CLOSE(6)
C====================	Read Flowrate data from FLUENT data ================C
c        OPEN(8,FILE='flowrate.dat',STATUS='UNKNOWN')
c	do 22 i=1,N
c	read(8,*) t,QR(i),QU(i)	  
c22	continue

	stop
	end
C-------------------------------------------------------------------------C
 	Recursive Subroutine impedance(s,Rroot,NA,NB,Z0,dummy)
	   USE MSFLIB
	COMMON/TEMPO/table(41,41)
	double complex ZD1,ZD2,ZD3,Z0,ZL,s,table,Zterm
	complex dummy
	real*8 pi
	Zterm=cmplx(0.,0.,8)
c	write(*,*) Zterm
	Rmin=0.003
	pi=datan(1.d0)*4.d0
	alpa=0.9
	beta=0.6
	R0=Rroot*alpa**(NA-1)*BETA**(NB-1)
c	write(5,*) NA,NB,R0
	if(R0.ge.Rmin) go to 10
            ZL=Zterm
	go to 100
c	write(*,*) table(NA,NB), dummy
10	   if(table(NA+1,NB).eq.dummy) go to 20
	      ZD1=table(NA+1,NB)
	   go to 25
20	      call impedance(s,Rroot,NA+1,NB,Z0,dummy)
c		  ZD1=Z0   
25	   continue
30	   if(table(NA,NB+1).eq.dummy) go to 40
	      ZD2=table(NA,NB+1)
	   go to 50
40	      call impedance(s,Rroot,NA,NB+1,Z0,dummy)
c		  ZD2=Z0
50	continue
	   ZL=ZD1*ZD2/(ZD1+ZD2)
100	CONTINUE
	call singlevesselimp(ZL,s,Z0,R0)
	write(5,300)
	write(5,*) ZD1,ZD2,ZL,Z0
300	format("ZD1=   ZD2=       ZL=      Z0=    ")
	table(NA,NB)=Z0
	return
	end
C===================================================C
	subroutine singlevesselimp(ZL,s,Z0,R0)
	   USE MSFLIB
	COMMON/TEMPO/table(41,41)
	Real*8 pi,del,C
	double complex Z0,ZL,s,ds,fun,table
	pi=datan(1.d0)*4.d0
	amu=0.0488
	gamma=2.
	alamda=20.   !p1311
	AL=alamda*R0
	rho=1.06
	C=1.5d0*pi*R0*R0/(2.E7*exp(-22.53*R0)+8.65E5)
	if(s.ne.0.) goto 10
	  Z0=ZL+2*(gamma+2)*amu*alamda/pi/r0**3
	goto 20
10	  del=2.d0*amu*(gamma+2.)/rho/r0/r0
	  ds=sqrt(pi*R0*R0/C/rho/s/(s+del))
	  call ctanh(AL/ds,fun)
	  z0=(ZL+fun/s/ds/C)/(s*ds*C*ZL*fun+1)
20	continue    
	write(5,30)
	write(5,*) s,C,del,ds,Z0
30	format("s=       C=       del=    ds=      Z0=    ")
      return
	end
C========================================================
	subroutine ctanh(x,y)
	double complex x,sum,y
	sum=exp(2*x)
	y=(sum-1)/(sum+1)
	return
	end