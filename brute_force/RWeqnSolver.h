#ifndef RWEQNSOLVER_H
#define RWEQNSOLVER_H

#include <complex>

using namespace std;

const complex<double> i(0.0, 1.0);

const int dim = 4;
const double l = 0.0;
const double omega = 0.4;
const double M = 0.5;

const double r0 = 1.001;
const double rInf = 100.0;

inline double Veff(double l, double M, double r)
{
	return (
			(1.0-2.0*M/r) * ( l*(l+1.0)/(r*r) + 2.0*M/(r*r*r) )
			);
}



#endif
