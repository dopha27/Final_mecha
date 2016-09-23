C	ImpedanceWeight	
      USE MSFLIB
	double complex Z(200),zn(200),z1,ZD1,ZD2,ZD3,Z0,zeta,theta,s
	double complex table,dummy
	COMMON/TEMPO/table(41,11)
	real*8 pi,aa,bb,eta,R0
	type (qwinfo) winfo
	Rroot=0.05
	dt=0.05
	N=20
	eps=1e-10
	eta=eps**(1/2/N)
	pi=datan(1.d0)*4.d0
	alpa=0.9
	beta=0.6
	dummy=(-1,0.)
c------------------------	Plot ----------------------
C        OPEN(5,FILE='cmplx1.dat',STATUS='UNKNOWN')
	open(2,file='user',title='Vorticity')
	istat=setwsizeqq(2,winfo)
	   call setviewport(0,0,400,300)
	   status = setwindow(.true.,-30,0,30,30)
  	    i1=SETCOLORRGB(255+255*256+256*256*255)   ! white Border
  	    dummy=RECTANGLE_w($gborder,-30,0,30,30)

	do 20 m=0,2*N-1
	do 10 i=1,40
	do 10 j=1,10
		table(i,j)=dummy
10	continue
	the=pi*m/N
	zeta=eta*cmplx(cos(the),sin(the),8)
c	write(*,*) zeta
	theta=zeta*zeta/2.-2*zeta+1.5
	s=theta/dt
 	call impedance(s,Rroot,1,1,Z0,dummy)
	Z(m+1)=Z0
20    continue
	PAUSE

        OPEN(6,FILE='cmplx.dat',STATUS='UNKNOWN')
	do 22 i=1,40
	do 22 j=1,10
	write(6,*) table(i,j)	  
22	continue
 
 	do 30 i=0,N
		z1=0.
	  do 25 j=0,2*N-1
		the=-i*j*pi/N
25	    z1=z1+Z(j+1)*cmplx(cos(the),sin(the),8)
	  zn(i+1)=z1*eta**(-i)/2/N
30	continue
      CLOSE(6)

	stop
	end
C-------------------------------------------------------------------------------------------------
 	Recursive Subroutine impedance(s,Rroot,NA,NB,Z0,dummy)
	   USE MSFLIB
	COMMON/TEMPO/table(41,11)
	double complex ZD1,ZD2,ZD3,Z0,ZL,s,table,dummy
   	real*8 wx,wy,u1,v1
 	type (wxycoord) wxy
	real*8 pi,R0
	Zterm=0.
	Rmin=0.003
	pi=datan(1.d0)*4.d0
	alpa=0.9
	beta=0.6
	R0=Rroot*alpa**(NA-1)*BETA**(NB-1)
C	write(*,*) NA,NB,R0
	if(R0.ge.Rmin) go to 10
            ZL=Zterm
	go to 100
c	write(*,*) table(NA,NB), dummy
10	   if(table(NA+1,NB).eq.dummy) go to 20
	      ZD1=table(NA+1,NB)
	   go to 25
20	      call impedance(s,Rroot,NA+1,NB,Z0,dummy)   
25	   continue
30	   if(table(NA,NB+1).eq.dummy) go to 40
	      ZD2=table(NA,NB+1)
	   go to 50
40	      call impedance(s,Rroot,NA,NB+1,Z0,dummy)
50	   ZL=ZD1*ZD2/(ZD1+ZD2)
100	CONTINUE
	call singlevesselimp(ZL,s,Z0,R0)
	table(NA,NB)=Z0
	do 101 NB=2,10
	do 101 NA=2,40
		wx=-(NA-1+NB-1)*dcos(25.d0*pi/180.)
		wy=NA-1
	if(table(NA,NB).eq.dummy) go to 101
		u1=-dcos(25.d0*pi/180.)
		v1=1
C  	    i1=SETCOLORRGB(i11)   ! white Border
 	    call moveto_w(wx,wy,wxy)
 	   status= lineto_w(wx+u1,wy+v1)
101	CONTINUE
	return
	end
C===================================================C
	subroutine singlevesselimp(ZL,s,Z0,R0)
	   USE MSFLIB
	COMMON/TEMPO/table(41,11),pi
	Real*8 pi,del,R0
	double complex Z0,ZL,s,ds,fun,table
	pi=datan(1.d0)*4.d0
	amu=0.0488
	gamma=2.
	alamda=20.   !p1311
	L=alamda*R0
	rho=1.06
	C=1.5*pi*R0*R0/(2.E7*exp(-22.53*R0)+8.65E5)
	if(s.ne.0.) goto 10
	  Z0=ZL+2*(gamma+2)*amu*alamda/pi/r0**3
	goto 20
10	  del=2*amu*(gamma+2)/rho/r0/r0
	  ds=sqrt(pi*R0*R0/C/rho/s/(s+del))
	  call ctanh(L/s,fun)
	  z0=(zl+fun/s/ds/C)/(s*ds*C*ZL*fun+1)
20	continue    
      return
	end
C========================================================
	subroutine ctanh(x,y)
	double complex x,sum,y
	sum=exp(2*x)
	y=(sum-1)/(sum+1)
	return
	end