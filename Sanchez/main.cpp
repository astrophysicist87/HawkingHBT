#include <iostream>
#include <complex>

#include "SanchezAsymptotics.h"

using namespace std;

int main(void)
{
	complex<double> phi_at_x = 0.0, phi_prime_at_x = 0.0;
	if (mode == 0)
		arbitrary_expansion();
	else
		horizon_expansion(phi_at_x, phi_prime_at_x);

	cout << phi_at_x.real() << "   " << phi_at_x.imag() << "   "
			<< phi_prime_at_x.real() << "   " << phi_prime_at_x.imag() << endl;

	return (0);
}
