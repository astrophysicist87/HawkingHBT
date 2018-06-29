#ifndef SANCHEZASYMPTOTICS_H
#define SANCHEZASYMPTOTICS_H

#include <iostream>
#include <complex>

using namespace std;

const complex<double> i(0.0, 1.0);

const int mode = 1;

void arbitrary_expansion();
void horizon_expansion(complex<double> & phi_at_x, complex<double> & phi_prime_at_x);

inline double disc(double t1, double t2)
{
	return ( 2.0*abs(t1-t2)/(t1+t2) );
}

#endif
