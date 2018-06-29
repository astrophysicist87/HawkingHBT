#include <iostream>
#include <complex>

using namespace std;

const complex<double> i(0.0, 1.0);

int main(void)
{
	double xs = 1.0;
	double b = 1.5;
	double y = b - xs;
	double x = 1.8;	//say

	const double l = 0;
	const int nmax = 100;

	//these are supposed to be arbitrary initial conditions
	complex<double> icC0 = 1.0, icC1 = -0.5;
	
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

	return (0);
}
