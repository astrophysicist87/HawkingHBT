#include <cstdlib>
#include <iostream>
#include <complex>

#include "SanchezAsymptotics.h"

using namespace std;

const bool first_coeff_only = true;

////////////////////////////////////////////////////////////////////////////////////////
void horizon_expansion( double x0,
						complex<double> & phi_at_x0,
						complex<double> & phi_prime_at_x0 )
{
	//double x0 = xs + stepsize;	//take initial step and get expansion here

	//set initial conditions at horizon
	complex<double> icd0 = 1.0;
	complex<double> icd1 = 0.0;
	
	//start off at n=2 to use Sanchez' Eq.(10)
	complex<double> dn = 0.0, dnm1 = 0.0, dnm2 = 0.0, dnm3 = 0.0;
	dnm1 = ( first_coeff_only ) ? icd0 : icd1;
	dnm2 = ( first_coeff_only ) ? 0.0 : icd0;

	complex<double> sum = icd0 + icd1*(x0-xs);	//incorporate first two terms immediately
	complex<double> sum2 = icd1;				//need this term for derivative
	double factor = ( first_coeff_only ) ? 1.0 : x0-xs;

	complex<double> prevsum = sum, prevsum2 = sum2;
	int in = ( first_coeff_only ) ? 1 : 2;
	do
	{
		double n = in;
		prevsum = sum;
		prevsum2 = sum2;

		complex<double> coeffn = (n-2.0*i*xs)*n*xs;
		complex<double> coeffnm1 = (n+l)*(n-l-1.0)
									+ 2.0*xs*xs
									- (2.0*n-1.0)*i*xs;
		complex<double> coeffnm2 = 3.0*xs;
		complex<double> coeffnm3 = 1.0;

		dn = - ( coeffnm1*dnm1
					+ coeffnm2*dnm2
					+ coeffnm3*dnm3 )
				/ coeffn;

/*
cout << "n = " << n << endl;
cout << "coeffnp2 = " << 0.0 << endl;
cout << "coeffnp1 = " << 0.0 << endl;
cout << "coeffn = " << coeffn << endl;
cout << "coeffnm1 = " << coeffnm1 << endl;
cout << "coeffnm2 = " << coeffnm2 << endl;
cout << "coeffnm3 = " << coeffnm3 << endl;
cout << "dnp2 = " << 0.0 << endl;
cout << "dnp1 = " << 0.0 << endl;
cout << "dn = " << dn << endl;
cout << "dnm1 = " << dnm1 << endl;
cout << "dnm2 = " << dnm2 << endl;
cout << "dnm3 = " << dnm3 << endl;
*/

		//sum2 defined with different factor!!!
		sum2 += dn * n * factor;	//term: d_n * n * (x-xs)^(n-1)
		factor *= x0-xs;
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
	complex<double> prefactor = exp(/*-i*xs*/ -i * xs * log(abs(x0-xs)));

	phi_at_x0 = prefactor * sum;
	phi_prime_at_x0 = prefactor * ( -i*xs*sum/(x0-xs) + sum2 );

	return;
}

////////////////////////////////////////////////////////////////////////////////////////
void arbitrary_expansion( double b, double x1,
							complex<double> phi_at_x0,
							complex<double> phi_prime_at_x0,
							complex<double> & phi_at_x1,
							complex<double> & phi_prime_at_x1 )
{
	double y = b - xs;
	if (y <= 0.0)
	{
		cerr << "y <= 0.0!!!  Something went wrong..." << endl;
		exit(8);
	}

	//set initial conditions at x==b in terms of phi(b) and phi'(b)
	complex<double> icC0 = exp(i*xs*log(y))*phi_at_x0;
	complex<double> icC1 = exp((i*xs-1.0)*log(y))
							* ( i*xs*phi_at_x0 + y * phi_prime_at_x0 );

	complex<double> Cn = icC0, Cnp1 = icC1, Cnp2 = 0.0;
	complex<double> Cnm1 = 0.0, Cnm2 = 0.0, Cnm3 = 0.0;

	complex<double> sum = icC0 + icC1*(x1-b);
	complex<double> sum2 = icC1;
	double factor = x1-b;

	int in = 0;
	complex<double> prevsum = sum, prevsum2 = sum2;
	do
	{
		double n = in;
		prevsum = sum;
		prevsum2 = sum2;

		complex<double> coeffnp2 = b*y*y*(n+1.0)*(n+2.0);
		complex<double> coeffnp1 = y*(n+1.0)*( y*(n+1.0) + b*(2.0*n+1.0-2.0*i*xs) );
		complex<double> coeffn = (b+2.0*y)*n*n
									+ (y-2.0*i*xs*(y+b))*n
									+ b*(b*b-xs*xs)
									- y*(l+i*xs);
		complex<double> coeffnm1 = (n-1.0)*(n-2.0*i*xs)+3.0*b*b - l/**(l+1.0)*/ - xs*(i+xs);
		complex<double> coeffnm2 = 3.0*b;
		complex<double> coeffnm3 = 1.0;
/*
cout << "n = " << n << endl;
cout << "coeffnp2 = " << coeffnp2 << endl;
cout << "coeffnp1 = " << coeffnp1 << endl;
cout << "coeffn = " << coeffn << endl;
cout << "coeffnm1 = " << coeffnm1 << endl;
cout << "coeffnm2 = " << coeffnm2 << endl;
cout << "coeffnm3 = " << coeffnm3 << endl;
cout << "Cnp2 = " << Cnp2 << endl;
cout << "Cnp1 = " << Cnp1 << endl;
cout << "Cn = " << Cn << endl;
cout << "Cnm1 = " << Cnm1 << endl;
cout << "Cnm2 = " << Cnm2 << endl;
cout << "Cnm3 = " << Cnm3 << endl;
*/

		Cnp2 = - ( coeffnp1*Cnp1
					+ coeffn*Cn
					+ coeffnm1*Cnm1
					+ coeffnm2*Cnm2
					+ coeffnm3*Cnm3 )
				/ coeffnp2;

		//sum2 defined with different factor!!!
		sum2 += Cnp2 * (n+2.0) * factor;
		factor *= x1-b;
		sum += Cnp2 * factor;

		Cnm3 = Cnm2;
		Cnm2 = Cnm1;
		Cnm1 = Cn;
		Cn = Cnp1;
		Cnp1 = Cnp2;

		//cout << n+2 << "   " << Cnp2 << "   " << Cnp1 / Cn << endl;
		//cout << n+2 << "   " << sum << endl;
		in++;
	} while (
				disc(prevsum.real(), sum.real()) > tolerance
				or disc(prevsum.imag(), sum.imag()) > tolerance
				or disc(prevsum2.real(), sum2.real()) > tolerance
				or disc(prevsum2.imag(), sum2.imag()) > tolerance
			);

	complex<double> prefactor = exp(-i * xs * log(abs(x1-xs)));

	phi_at_x1 = prefactor * sum;
	phi_prime_at_x1 = prefactor * ( -i*xs*sum/(x1-xs) + sum2 );

	return;
}


