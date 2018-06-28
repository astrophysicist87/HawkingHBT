#ifndef YRT_H
#define YRT_H

#include <stdio.h>
#include <iostream>
#include <complex>
#include <vector>

#include <gsl/gsl_sf_gamma.h>

using namespace std;

namespace Yrt_NS
{
	const complex<double> i(0.0,1.0);

	extern int nmax, rmax, smax, l, n_x_pts;

	void initialize_A3D_starting_values();
	void set_B1D(int n);
	void set_C1D(int n);
	void set_A3D(int n, int r, int s);
	void set_B3D(int n, int t, int s);
	void set_C3D(int n, int t, int s);
	void get_Yrt(
			vector<complex<double> > * Yrt_in,
			vector<complex<double> > * ddx_Yrt_in,
			vector<double> * x_pts_in,
			int nmax_in, int rmax_in, int l_in);

	inline double pp(int r0)
	{
		if (r0 > l || r0 < 0)
			return (0.0);
		else
		{
			int r = l - r0;	//corrects for typo in Persides' paper
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
		else if (r > 95)
		{
			double prefactor = pow(-1.0, l+1) / 2.0;
			double ln_r_fact = gsl_sf_lnfact(r);
			double ln_twolprp1_fact = gsl_sf_lnfact(2*l+r+1);
			double ln_lpr_fact = gsl_sf_lnfact(l+r);
			double ln_result_no_pref = 2.0*ln_lpr_fact - ln_r_fact - ln_twolprp1_fact;

			return ( prefactor * exp( ln_result_no_pref ) );
		}
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
		if (s > 165 or t > 165)
		{
			double prefactor = pow(-1.0, s-t) / pow(r+1.0, s-t+1.0);
			double ln_ratio = gsl_sf_lnfact(s) - gsl_sf_lnfact(t);
			return (prefactor*exp(ln_ratio));
		}
		else if ( r != -1 and 0 <= t and t < s+1 )
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

	/*inline int indexer_A4D(int n, int r, int s, int t)
	{
		//range of t is 0<=t<=s ==> just use smax
		return ( ( ( n * rmax + r ) * smax + s ) * smax + t );
	}*/

	inline int indexer_Yrt(int r, int t, int ix)
	{
		//range of t is 0<=t<=s ==> just use smax
		return ( ( r * smax + t ) * n_x_pts + ix );
	}

}

#endif
