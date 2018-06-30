#include <cstdlib>
#include <iostream>
#include <complex>
#include <fenv.h>

#include "SanchezAsymptotics.h"

using namespace std;

int main(void)
{
	feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

	complex<double> phi_at_x0 = 0.0, phi_prime_at_x0 = 0.0;

	double x0 = xs;	//take initial step and get expansion here
	const int ninit = 1;
	//for (int tmp = 0; tmp < 100000-1; ++tmp)
	for (int tmp = 0; tmp < ninit; ++tmp)
	//for (int tmp = 0; tmp < 1; ++tmp)
	{
		x0 += stepsize;

		horizon_expansion(x0, phi_at_x0, phi_prime_at_x0);

		cout << x0 << "   "
				<< phi_at_x0.real() << "   "
				<< phi_at_x0.imag() << "   "
				<< phi_prime_at_x0.real() << "   "
				<< phi_prime_at_x0.imag() << endl;
	}
if (ninit>1) return (0);
	int nsteps = ninit+1;
	double xold = x0, xnew = x0 + stepsize;
	for (int step = 0; step < nsteps; ++step)
	{
		complex<double> phi_at_x = phi_at_x0;
		complex<double> phi_prime_at_x = phi_prime_at_x0;

		arbitrary_expansion(xold, xnew,
							phi_at_x0, phi_prime_at_x0,
							phi_at_x, phi_prime_at_x);

		cout /*<< "(arb): " */<< xnew << "   "
				<< phi_at_x.real() << "   "
				<< phi_at_x.imag() << "   "
				<< phi_prime_at_x.real() << "   "
				<< phi_prime_at_x.imag() << endl;

		//update to next expansion point
		phi_at_x0 = phi_at_x;
		phi_prime_at_x0 = phi_prime_at_x;
		xold += stepsize;
		xnew += stepsize;
	}

	return (0);
}
