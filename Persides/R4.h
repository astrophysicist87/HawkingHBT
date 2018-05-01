#include <stdio.h>
#include <iostream>
#include <complex>
#include <vector>

#include <gsl/gsl_sf_gamma.h>

using namespace std;

const complex<double> i(0.0,1.0);

const int l = 0;

inline double pp(int r)
{
	if (r > l || r < 0)
		return (0.0);
	else
	{
		double r_fact = gsl_sf_fact(r);
		double lpr_fact = gsl_sf_fact(r);
		double lmr_fact = gsl_sf_fact(r);

		return ( pow(-1.0, r)*lpr_fact / ( r_fact*r_fact*lmr_fact ) );
	}
}

inline double qp(int r)
{
	if (r < 0)
		return (0.0);
	else
	{
		double twor_fact = 2.0*gsl_sf_fact(r);
		double twolprp1_fact = gsl_sf_fact(r);
		double lpr_fact = gsl_sf_fact(r);

		return ( pow(-1.0, l+1)*lpr_fact*lpr_fact / ( twor_fact*twolprp1_fact ) );
	}
}

inline int indexer_A3D(int n, int r, int s)
{
		return ( ( n * rmax + r ) * smax + s );
}
