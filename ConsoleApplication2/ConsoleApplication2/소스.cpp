#include <complex>
#include <iostream>
#include <math.h>
#include <vector>

/*
int main()
{
	using namespace std;

	complex <double> c1(1.0, 1.0);
	complex <double> c2(4.0, 3.0);
	cout << "The complex number c1 = " << c1 << endl;

	double dr1 = c1.real();
	cout << "The real part of c1 is c1.real ( ) = "
		<< dr1 << "." << endl;

	double di1 = c1.imag();
	cout << "The imaginary part of c1 is c1.imag ( ) = "
		<< di1 << "." << endl;

	complex <double>z1 ;
	z1 = c1 * c2;
	cout << "The complex number z1 = " << z1 << endl;

}
*/

#define PI 3.141592654   // #define I _Complex_I #undef I

double rho = 1.06; /*rho=density, gamma=velocity profile, mu=blood viscosity*/
double mu = 0.0488;
double gamma = 2;

double SINGLEVESSELIMP_imag;



using namespace std;

double IMPEDANCE(double r_root, complex<double> s, int Na, int Nb, complex <double> table[38][38] , complex<double> Z[200] , int m);
double SINGLEVESSELIMP(complex<double> ZL, complex <double> s, double r0, complex<double> *Z0);

double *IMPEDANCEWEIGHTS(double r_root, double dt, double eps)
{
	//_C_double_complex zta, Sigma, a, two;

	double x;
	complex <double> table[38][38];
	int n, m;
	int k, j;
	int N = ceil(1 / dt); 	   // ceil is smallest following integer to 1/	
	
	complex <double> Zn_final[100] = {0.0};
	complex <double> x1 = 0.0;
	double eta = pow(eps, (1 / (2 * N)));	// check again
	complex <double> Z[200] = { 0.0 };

	for (m = 0; m < (2 * N - 1); m++)
	{
		complex <double> a (0.0, PI*m / N);
		complex <double> zta;
		zta = eta*exp(a);
		complex <double> two (2,0);
		complex <double> Sigma;
		Sigma = zta*zta/2.0 - 2.0*zta+1.5;

		x = sqrt(-1.0);
		for (k = 0; k<38; k++) {
			for (j = 0; j<38; j++) {
				table[k][j] = x;
			}
		}

		IMPEDANCE(r_root, Sigma / dt, 1, 1, table , Z , m);
		//	 printf("%lf + %lf i \n",a1,a2);
		

	
	}
	
	for (n = 0; n < N; n++)	{
		for (m = 0; m < (2 * N - 1); m++) {
			x1 += Z[m]*exp();
		
		}
		x1 = x1*pow(eta, -n) / (2 * N);
		

		//printf("Zn = %lf + %lf i \n", x1, x2);

	}
	
	return 0;
}


double IMPEDANCE(double r_root, complex<double> s, int Na, int Nb, complex <double> table[38][38], complex<double> Z[200], int m)
{
	/* F recursive function */
	/* Output: Laplace transform of tree impedance */

	double alpha = 0.91;
	double beta = 0.58;
	double r_min = 0.003;
	double Z_term = 0;	// ???
	complex <double> ZL, ZD1, ZD2;
	complex <double> Z0;

	double r0 = r_root*pow(alpha, (Na - 1))*pow(beta, (Nb - 1));

	if (r0 < r_min)
	{
		ZL = Z_term;
	}
	else
	{
		if (imag(table[Na + 1][Nb]) != 0)
		{
			ZD1 = IMPEDANCE(r_root, s, Na + 1, Nb, table);
		}
		else
		{
			ZD1 = table[Na + 1][Nb];
		}

		if (imag(table[Na][Nb + 1]) != 0)
		{
			ZD2 = IMPEDANCE(r_root, s, Na, Nb + 1, table);
		}
		else
		{
			ZD2 = table[Na][Nb + 1];
		}

		ZL = ZD1 * ZD2 / (ZD1 + ZD2);
		//	printf("ZD1=%f ||| ZD2=%f ||| ZL=%f \n", ZD1, ZD2, ZL);
	}

	// Z0를 a1과 a2로 나타내서 만들기.

	 
	SINGLEVESSELIMP(ZL, s, r0, &Z0);
	//printf("Na = %d, Nb = %d, Z0 = %f\n", Na, Nb, Z0);
	//printf("%g", Z0);
	table[Na][Nb] = Z0;

	Z[m] = Z0;
																					// Z 0 -> return 해주기
	return real(Z0);
}


double SINGLEVESSELIMP(complex<double> ZL, complex <double> s, double r0, complex<double> *Z0)
{
	/* Output: Z0, Laplace transform of the impedance at x = 0 given ZL, its value at x = L */

	complex <double> ds;
	double del;
	double lamda, L, C, Eh, k1, k2, k3;/*lamda=length/radius ratio, C=compliance, Eh=elasticity*/
	double A1, B1, C1, D1;
	double e, COS, SIN;

	if (r0 > 0.000250)
	{
		lamda = 10;
		L = lamda*r0;
	}
	if ((r0 >= 0.000050) && (r0 <= 0.000250))
	{
		lamda = 20;
		L = lamda*r0;
	}
	if (r0 < 0.000050)
	{
		lamda = 30;
		L = lamda*r0;
	}
	k1 = 2 * pow(10.0, 7);
	k2 = -22.53;
	k3 = 8.65*pow(10.0, 5);
	Eh = k1*exp(k2*r0) + k3;
	C = 1.5*PI*r0*r0*r0 / Eh;


	if (real(s)== 0 && imag (s)== 0)
	{
		*Z0 = ZL + 2 * (gamma + 2)*mu*lamda / (PI*r0*r0*r0);
		
	}
	else
	{
		del = 2 * mu*(gamma + 2) / (rho*r0*r0);
		ds = sqrt((PI*r0*r0) / (C*rho*s*(s + del)));
		*Z0 = ((ZL + tanh(L/ds)/s*ds*C)/(s*ds*C*ZL*tanh(L/ds)+1.0));

	}
	cout << ds << endl;
	//printf("%lf  +  %lf i \n", real(ds), imag(ds));
	return real(ZL);
}

void main()
{
	double dt = 0.05;
	double r_root = 0.4;
	double eps = pow(10.0, -10);
	double*Z;


	Z = IMPEDANCEWEIGHTS(r_root, dt, eps);

}
