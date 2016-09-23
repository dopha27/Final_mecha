C	Main Program for Impedance BC (Pressure)
      USE MSFLIB
	COMMON/TEMPO/table(41,41)
	double complex zn(200),table
c	complex table
	dimension QR(200),QU(200),PBC(200)
        OPEN(3,FILE='flowrate.dat',STATUS='UNKNOWN')
	Pterm=45.
	Rroot=0.05
	dt=0.05
	N=20
	pause
	call ImpedanceWeight(zn,Rroot,dt,N)	
	do 10 i=1,2*N+1
	read(3,*) a,b,c
	QR(i)=b
	QU(i)=c
10	continue
	do 30 i=1,N+19
	pibc=Pterm
	do 20 k=0,i
	pibc=pibc+Real(zn(k+1))*QR(i-k+1)
20	continue
30	PBC(i)=pibc
	write(*,40) (i,PBC(i), i=1,N+19)
40	format(2x,'i=  ',i2,'   P0=  ',f10.6)
	stop
	end
C------------------------------------------------------
	subroutine ImpedanceWeight(zn,Rroot,dt,N)
C--------------------------------------------------------		
      USE MSFLIB
	COMMON/TEMPO/table(41,41)
	double complex Z(200),zn(200),z1,ZD1,ZD2,ZD3,Z0,zeta,theta,s,table
	complex dummy
	real*8 pi,eta
	eps=1e-10
	eta=eps**(1.d0/2/N)
	pi=datan(1.d0)*4.d0
	alpa=0.9
	beta=0.6
	dummy=(-3.,0.)
        OPEN(5,FILE='cmplx1.dat',STATUS='UNKNOWN')

	do 20 m=0,2*N-1
      	do 10 i=1,40
	    do 10 j=1,40
		    table(i,j)=dummy
10	continue
	  the=pi*m/N
	  zeta=eta*cmplx(cos(the),sin(the),8)
	  theta=zeta*zeta/2.-2.*zeta+1.5
	  s=theta/dt
	write(5,*) s
 	  call impedance(s,Rroot,1,1,Z0,dummy)
	  Z(m+1)=Z0
20    continue
	close(5)
        OPEN(6,FILE='cmplx.dat',STATUS='UNKNOWN')
c	do 22 j=1,7
c	do 22 i=1,30
c	write(6,*) table(i,j)	  
c22	continue
 
 	do 30 i=0,N
		z1=(0.,0.)
	  do 25 j=0,2*N-1
		the=-i*j*pi/N
25	    z1=z1+Z(j+1)*cmplx(cos(the),sin(the),8)
	  zn(i+1)=z1*eta**(-i)/2/N
	write(6,*) i,zn(i+1)
30	continue
      CLOSE(6)
	return
	end
C-------------------------------------------------------------------------C
 	Recursive Subroutine impedance(s,Rroot,NA,NB,Z0,dummy)
	   USE MSFLIB
	COMMON/TEMPO/table(41,41)
	double complex ZD1,ZD2,ZD3,Z0,ZL,s,Zterm,table
	complex dummy
	real*8 pi
	Zterm=cmplx(0.,0.,8)
	Rmin=0.003
	pi=datan(1.d0)*4.d0
	alpa=0.9
	beta=0.6
	R0=Rroot*alpa**(NA-1)*BETA**(NB-1)
	if(R0.ge.Rmin) go to 10
            ZL=Zterm
	go to 100
10	   if(table(NA+1,NB).eq.dummy) go to 20
	      ZD1=table(NA+1,NB)
	   go to 30
20	      call impedance(s,Rroot,NA+1,NB,Z0,dummy)
		  ZD1=Z0   
30	   continue
	   if(table(NA,NB+1).eq.dummy) go to 40
	      ZD2=table(NA,NB+1)
	   go to 50
40	      call impedance(s,Rroot,NA,NB+1,Z0,dummy)
		  ZD2=Z0
50	continue
	   ZL=ZD1*ZD2/(ZD1+ZD2)
100	CONTINUE
	call singlevesselimp(ZL,s,Z0,R0)
c	write(5,300)
c	write(5,*) ZD1,ZD2,ZL,Z0
300	format("ZD1=   ZD2=       ZL=      Z0=    ")
	table(NA,NB)=Z0
	return
	end
C===================================================C
	subroutine singlevesselimp(ZL,s,Z0,R0)
	   USE MSFLIB
	Real*8 pi,del,C
	double complex Z0,ZL,s,ds,fun
	pi=datan(1.d0)*4.d0
	amu=0.0488
	gamma=2.
	alamda=30.   !p1311
	if(R0.ge.0.005) alamda=20.
	if(R0.ge.0.025) alamda=10.
	AL=alamda*R0
	rho=1.06
	C=1.5d0*pi*R0*R0/(2.E7*exp(-22.53*R0)+8.65E5)
	if(s.ne.0.) goto 10
	  Z0=ZL+2*(gamma+2)*amu*alamda/pi/r0**3
	goto 20
10	  del=2.d0*amu*(gamma+2.)/rho/r0/r0
	  ds=sqrt(pi*R0*R0/C/rho/s/(s+del))
	  call ctanh(AL/ds,fun)
c		fun=tanh(AL/ds)
	  z0=(ZL+fun/s/ds/C)/(s*ds*C*ZL*fun+1)
20	continue    
c	write(5,30)
c	write(5,*) s,C,del,ds,Z0
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