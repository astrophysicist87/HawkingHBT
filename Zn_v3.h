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

//some parameters
const double xmin = 1.0;
const double xmax = 10.0;
const int nmax = 2;
const int lmax = 1;
const int n_x_pts = 1001;
const int order = n_x_pts - 1;
const complex<double> i(0,1);

//function prototypes
inline int indexer(int il, int ix);
inline int indexer(int il, int in, int ix);

void get_hankel(int l, double x,
				complex<double> * h1l,
				complex<double> * h2l);
void convert_array_to_Chebyshev_coeffs(double a_in, double b_in,
										vector<double> values_on_nodes,
										gsl_cheb_series * coefficients);
void set_Z0(int l);
void set_ddx_Zn(int il, int in);
void set_ddx_x_ddx_Zn(int il, int in);
void set_integrand_1_coeffs(int il, int in,
							gsl_cheb_series * csIntegral_re,
							gsl_cheb_series * csIntegral_im);
void set_integrand_2_coeffs(int il, int in,
							gsl_cheb_series * csIntegral_re,
							gsl_cheb_series * csIntegral_im);
void set_Zn(int l, int in);

#endif
