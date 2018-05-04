#include <stdio.h>
#include <iostream>
#include <complex>
#include <vector>

#include <gsl/gsl_sf_gamma.h>

using namespace std;

vector<complex<double> > B1D, C1D;
vector<complex<double> > A3D, B3D, C3D;
vector<complex<double> > A4D;

const complex<double> i(0.0,1.0);
const int rmax = 80, nmax = 3, smax = 3;
const int tmax = rmax;
const int l = 2;

void initialize_A3D_starting_values();
void set_B1D(int n);
void set_C1D(int n);
void set_A3D(int n, int r, int s);
void set_B3D(int n, int t, int s);
void set_C3D(int n, int t, int s);

inline double pp(int r)
{
	if (r > l || r < 0)
		return (0.0);
	else
	{
		double r_fact = gsl_sf_fact(r);
		double lpr_fact = gsl_sf_fact(l+r);
		double lmr_fact = gsl_sf_fact(l-r);

		return ( pow(-1.0, r)*lpr_fact / ( r_fact*r_fact*lmr_fact ) );
	}
}

inline double qp(int r)
{
	if (r < 0)
		return (0.0);
	else
	{
		double r_fact = gsl_sf_fact(r);
		double twolprp1_fact = gsl_sf_fact(2*l+r+1);
		double lpr_fact = gsl_sf_fact(l+r);

		return ( pow(-1.0, l+1)*lpr_fact*lpr_fact / ( 2.0*r_fact*twolprp1_fact ) );
	}
}

inline double Al(int il)
{
	double l_fact = gsl_sf_fact(l);
	double twol_fact = gsl_sf_fact(2*l);
	double twolpone_fact = gsl_sf_fact(2*l+1);

	double num = pow(-2.0,l)*pow(l_fact, 3.0);
	double den = twol_fact*twolpone_fact;

	return ( num / den );
}

inline double D3D(int r, int s, int t)	//N.B. - not the same as r, s, t elsewhere!
{
	if ( r != -1 and 0 <= t and t < s+1 )
	{
		double num = pow(-1.0, s-t)*gsl_sf_fact(s);
		double den = pow(r+1.0, s-t+1.0)*gsl_sf_fact(t);
		return (num/den);
	}
	else if ( r == -1 and t == s+1 )
		return ( 1.0 / (s + 1.0) );
	else
		return (0.0);
}

inline int indexer_A3D(int n, int r, int s)
{
		return ( ( n * rmax + r ) * smax + s );
}

inline int indexer_B3D(int n, int t, int s)
{
		return ( ( n * rmax + t ) * smax + s );
}

inline int indexer_C3D(int n, int t, int s)
{
		return ( ( n * rmax + t ) * smax + s );
}



