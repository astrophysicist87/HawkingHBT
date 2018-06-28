#ifndef XN_H
#define XN_H

#include <stdio.h>
#include <iostream>
#include <complex>
#include <vector>

#include <gsl/gsl_sf_gamma.h>

using namespace std;

namespace Xn_NS
{
	const complex<double> i(0.0,1.0);

	extern int nmax, l, n_x_pts;

	void get_Xn(
			vector<complex<double> > * Xn_in,
			vector<double> * x_pts_in,
			int nmax_in, int l_in);

	void set_Pl_arg_and_Ql_arg(int l_local);
	void set_Xt0_and_Xt1(double sigma_l);

	inline int indexer_Xn(int n, int ix)
	{
		return ( n * n_x_pts + ix );
	}

	inline int indexer_Xtn(int n, int ix)
	{
		return ( n * n_x_pts + ix );
	}

	inline double pp(int k, int l_local)
	{
		if (k > l_local || k < 0)
			return (0.0);
		else
		{
			//int r = l_local - k;	//corrects for typo in Persides' paper
			int r = k;	//corrects for typo in Persides' paper
			double r_fact = gsl_sf_fact(r);
			double lpr_fact = gsl_sf_fact(l_local+r);
			double lmr_fact = gsl_sf_fact(l_local-r);

			return ( pow(-1.0, r)*lpr_fact / ( r_fact*r_fact*lmr_fact ) );
		}
	}

	inline double eval_poly(double coeffs [], int order, double x)
	{
		double factor = 1.0, result = 0.0; 
		for(int term = 0; term <= order; term++) {
			result += coeffs[term] * factor;
			factor *= x;
		}
		return result;
	}
}

#endif
