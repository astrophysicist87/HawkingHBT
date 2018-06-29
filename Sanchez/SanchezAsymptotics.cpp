#include <iostream>
#include <complex>

#include "SanchezAsymptotics.h"

using namespace std;

void arbitrary_expansion()
{
	double xs = 1.0;
	double b = 1.5*xs;
	double y = b - xs;
	double x = 2.0*xs;	//say

	const double l = 0;
	const int nmax = 1000;

	//these are supposed to be arbitrary initial conditions
	complex<double> icC0 = 1.0;
	complex<double> icC1 = i*xs*icC0/(b-xs);	//sets derivative to zero at horizon
	
	complex<double> Cn = icC0, Cnp1 = icC1, Cnp2 = 0.0;
	complex<double> Cnm1 = 0.0, Cnm2 = 0.0, Cnm3 = 0.0;

	complex<double> sum = icC0 + icC1*(x-b);
	double factor = x-b;

	for (int in = 0; in < nmax; ++in)
	{
		double n = in;
		factor *= x-b;
		complex<double> coeffnp2 = b*y*y*(n+1.0)*(n+2.0);
		complex<double> coeffnp1 = y*(n+1.0)*( y*(n+1.0) + b*(2.0*n+1.0-2.0*i*xs) );
		complex<double> coeffn = (b+2.0*y)*n*n
									+ (y-2.0*i*xs*(y+b))*n
									+ b*(b*b-xs*xs)
									- y*(l+i*xs);
		complex<double> coeffnm1 = (n-1.0)*(n-2.0*i*xs)+3.0*b*b - l - xs*(i+xs);
		complex<double> coeffnm2 = 3.0*b;
		complex<double> coeffnm3 = 1.0;

		Cnp2 = - ( coeffnp1*Cnp1
					+ coeffn*Cn
					+ coeffnm1*Cnm1
					+ coeffnm2*Cnm2
					+ coeffnm3*Cnm3 )
				/ coeffnp2;

		sum += Cnp2 * factor;

		Cnm3 = Cnm2;
		Cnm2 = Cnm1;
		Cnm1 = Cn;
		Cn = Cnp1;
		Cnp1 = Cnp2;

		//cout << n+2 << "   " << Cnp2 << "   " << Cnp1 / Cn << endl;
		cout << n+2 << "   " << sum << endl;
	}
	return;
}



void horizon_expansion(complex<double> & phi_at_x, complex<double> & phi_prime_at_x)
{
	const double tolerance = 1.e-10;

	double xs = 1.0;
	double x = 1.5*xs;	//say

	const double l = 0;
	const int nmax = 25;

	//set initial conditions
	complex<double> icd0 = 1.0;
	complex<double> icd1 = 0.0;
	
	//start off at n=2 to use Sanchez' Eq.(10)
	complex<double> dn = 0.0, dnm1 = icd1, dnm2 = icd0, dnm3 = 0.0;

	complex<double> sum = icd0 + icd1*(x-xs);	//incorporate first two terms immediately
	complex<double> sum2 = icd1;				//need this term for derivative
	double factor = x-xs;

	complex<double> prevsum = sum, prevsum2 = sum2;
	int in = 2;
	//for (int in = 2; in < nmax; ++in)
	do
	{
		double n = in;
		prevsum = sum;
		prevsum2 = sum2;

		complex<double> coeffn = (n-2.0*i*xs)*n*xs;
		complex<double> coeffnm1 = (n+l)*(n-l-1.0) + 2.0*xs*xs - (2.0*n-1.0)*i*xs;
		complex<double> coeffnm2 = 3.0*xs;
		complex<double> coeffnm3 = 1.0;

		dn = - ( coeffnm1*dnm1
					+ coeffnm2*dnm2
					+ coeffnm3*dnm3 )
				/ coeffn;

		//sum2 defined with different factor!!!
		sum2 += dn * n * factor;	//term: d_n * n * (x-xs)^(n-1)
		factor *= x-xs;
		sum += dn * factor;			//term: d_n * (x-xs)^n

		dnm3 = dnm2;
		dnm2 = dnm1;
		dnm1 = dn;

		//cout << n << "   " << sum.real() << "   " << sum.imag() << "   "
		//		<< sum2.real() << "   " << sum2.imag() << endl;
		in++;
	} while (
				disc(prevsum.real(), sum.real()) > tolerance
				or disc(prevsum.imag(), sum.imag()) > tolerance
				or disc(prevsum2.real(), sum2.real()) > tolerance
				or disc(prevsum2.imag(), sum2.imag()) > tolerance
			);

	//cout << "FINAL:   " << sum.real() << "   " << sum.imag() << "   "
	//		<< sum2.real() << "   " << sum2.imag() << endl;
	complex<double> prefactor = exp(-i * xs * log(abs(x-xs)));

	phi_at_x = prefactor * sum;
	phi_prime_at_x = prefactor * ( -i*xs*sum/(x-xs) + sum2 );

	return;
}
