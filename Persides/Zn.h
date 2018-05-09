#ifndef ZN_H
#define ZN_H

#include <stdio.h>
#include <iostream>
#include <complex>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_sf_bessel.h>

#include "chebyshev_library.h"
#include "gauss_quadrature.h"

namespace Zn_NS
{
	extern int n_x_pts, nmax, lmax;

	//some parameters
	const complex<double> i(0,1);

	//function prototypes
	void linspace(vector<double> & x, double a, double b);

	void get_hankel(int l, double x,
					complex<double> * h1l,
					complex<double> * h2l);
	void convert_array_to_Chebyshev_coeffs(double a_in, double b_in,
											vector<double> values_on_nodes,
											gsl_cheb_series * coefficients);
	void set_Z0(int l);
	void set_ddx_Zn(int il, int in);
	void set_d2dx2_Zn(int il, int in);
	void set_ddx_x_ddx_Zn(int il, int in);
	void set_integrand_1_coeffs(int il, int in,
								gsl_cheb_series * csIntegral_re,
								gsl_cheb_series * csIntegral_im);
	void set_integrand_2_coeffs(int il, int in,
								gsl_cheb_series * csIntegral_re,
								gsl_cheb_series * csIntegral_im);
	void set_Zn(int l, int in);

	void compute_running_integration(
			vector<complex<double> > * fpts,
			vector<complex<double> > * integ_fpts);
	void get_all_Zn_and_ddx_Zn(
			vector<complex<double> > * Zn_array,
			vector<complex<double> > * ddx_Zn_array,
			vector<double> * xpts,
			int nmax, int lmax);

	//header file definitions
	inline int indexer(int il, int ix)
	{
		return( il * n_x_pts + ix );
	}

	inline int indexer(int il, int in, int ix)
	{
		return( ( il * nmax + in ) * n_x_pts + ix );
	}
}

#endif
