#ifndef RWEQNSOLVER_H
#define RWEQNSOLVER_H

#include <complex>

using namespace std;

const complex<double> i(0.0, 1.0);

const int dim = 4;
const double l = 0.0;
const double omega = 0.4;
const double M = 0.5;

const double r0 = 1.0+1.e-5;
const double rInf = 20.0;

const double relative_precision = 1.e-12;

inline double Veff(double l, double M, double r)
{
	return (
			(1.0-2.0*M/r) * ( l*(l+1.0)/(r*r) + 2.0*M/(r*r*r) )
			);
}


inline void solve_2D_simult_eqns(
				complex<double> Rin_at_rInf, complex<double> Rin_prime_at_rInf,
				complex<double> Rout_at_rInf, complex<double> Rout_prime_at_rInf,
				complex<double> Rc_at_rInf, complex<double> Rc_prime_at_rInf,
				complex<double> & Ain, complex<double> & Aout )
{
	complex<double> det = Rin_at_rInf * Rout_prime_at_rInf - Rout_at_rInf * Rin_prime_at_rInf;
	Ain = ( Rc_at_rInf * Rout_prime_at_rInf - Rout_at_rInf * Rc_prime_at_rInf ) / det;
	Aout = - ( Rc_at_rInf * Rin_prime_at_rInf - Rin_at_rInf * Rc_prime_at_rInf ) / det;
	
	return;
}



#endif
