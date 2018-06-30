#ifndef SANCHEZASYMPTOTICS_H
#define SANCHEZASYMPTOTICS_H

#include <cstdlib>
#include <iostream>
#include <complex>

using namespace std;

const complex<double> i(0.0, 1.0);

const double stepsize = 0.00001;
//const double stepsize = 0.1;
const double tolerance = 1.e-10;

const double xs = 1.0;
const double l = 0;

void arbitrary_expansion( double b, double x1,
							complex<double> phi_at_x0,
							complex<double> phi_prime_at_x0,
							complex<double> & phi_at_x1,
							complex<double> & phi_prime_at_x1 );
void horizon_expansion( double x0,
						complex<double> & phi_at_x0,
						complex<double> & phi_prime_at_x0 );

inline double disc(double t1, double t2)
{
	return ( 2.0*abs(t1-t2)/(t1+t2) );
}

/*USAGE:
debugger(__LINE__, __FILE__);
*/
void inline debugger(int cln, const char* cfn)
{
	cerr << "You made it to " << cfn << ":" << cln << "!" << endl;
	return;
}


#endif
