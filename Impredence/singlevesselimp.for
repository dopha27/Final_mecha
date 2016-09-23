      USE MSFLIB    
c      USE complex.h 
	Real*8 pi,del
	double complex Z0,ZL,s,ds,fun
c      int main() {
	pi=datan(1.d0)*4.d0
	amu=0.0488
	gamma=2.
	alamda=30.
	L=6.
	s=(2.,3.)
	R0=0.01
	rho=1.06
	ZL=(20.,10.)
	C=1.5*pi*R0*R0/(2.E7*exp(-22.53*R0)+8.65E5)
	if(s.eq.0.) goto 10
	  Z0=ZL+2*(gamma+2)*amu*alamda/pi/r0**3
	goto 20
10	  del=2*amu*(gamma+2)/rho/r0/r0
	  ds=sqrt(pi*r0*r0/C/rho/s/(s+del))
	  call ctanh(L/s,fun)
	  z0=(zl+fun/s/ds/C)/(s*ds*C*ZL*fun+1)
    
20      write(*,*) Z0

c      double complex quotient = z1 / z2;
c      printf("The quotient: Z1 / Z2 = %.2f %+.2fi\n", creal(quotient), cimag(quotient));

c      double complex conjugate = conj(z1);
c      printf("The conjugate of Z1 = %.2f %+.2fi\n", creal(conjugate), cimag(conjugate));

      stop
	end

	subroutine ctanh(x,y)
	double complex x,sum,y
	sum=exp(2*x)
	y=(sum-1)/(sum+1)
	return
	end