#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_sf_bessel.h>

#include "chebyshev_library.h"
#include "gauss_quadrature.h"
#include "Zn.h"

using namespace std;
using namespace csf;

namespace Zn_NS
{
	int nmax, lmax, n_x_pts;

	vector<double> x_pts;
	vector<complex<double> > Zn_array;
	vector<complex<double> > ddx_Zn_array;
	vector<complex<double> > d2dx2_Zn_array;
	vector<complex<double> > ddx_x_ddx_Zn_array;

	void linspace(vector<double> & x, double a, double b)
	{
	//returns vector x of n linearly spaced points, with a and b as endpoints
		int n = x.size();
		double Del_x = (b-a)/(double)(n-1);

		//assume x has length n already
		for (int i = 0; i < n; i++)
			x[i] = a + Del_x * (double)i;

		return;
	}

	void get_hankel(int l, double x, complex<double> * h1l, complex<double> * h2l)
	{
		//cout << "get_hankel(): " << l << "   " << x << endl;
		complex<double> jl = gsl_sf_bessel_jl(l, x);
		complex<double> yl = gsl_sf_bessel_yl(l, x);
		*h1l = jl + i * yl;
		*h2l = jl - i * yl;
		return;
	}


	void set_Z0(int l)
	{
		//set Z0
		for (int ix = 0; ix < n_x_pts; ++ix)
		{
			double x = x_pts[ix];
			complex<double> eix = exp(i*x), emix = exp(-i*x);
			complex<double> h1l, h2l, h1lp1, h2lp1;
			get_hankel(l, x, &h1l, &h2l);
			get_hankel(l+1.0, x, &h1lp1, &h2lp1);

			complex<double> coeff1(-2.0*x*x+l*l,x*(2.0*l+1.0));
			complex<double> coeff2 = x*(1.0-2.0*i*x);

			//set Z0
			Zn_array[indexer(l, 0, ix)] = emix * h1l;

			//be careful about which of these is which!!!
			//set Z0 derivative combinations in integral expressions
			//Z0_deriv1_array[indexer(l, ix)] = eix * ( coeff1 * h2l + coeff2 * h2lp1 ) / x;
			//Z0_deriv2_array[indexer(l, ix)] = eix * ( coeff1 * h1l + coeff2 * h1lp1 ) / x;
		}

		//set first derivative
		set_ddx_Zn(l, 0);

		//set second derivative
		set_d2dx2_Zn(l, 0);

		//set derivative combination
		set_ddx_x_ddx_Zn(l, 0);

		return;
	}

	void set_ddx_Zn(int il, int in)
	{
		double h = x_pts[1] - x_pts[0];

		//approximate first derivative at lower endpoint
		ddx_Zn_array[indexer(il, in, 0)]
			= ( - 3.0 * Zn_array[indexer(il, in, 0)]
				+ 4.0 * Zn_array[indexer(il, in, 1)]
				- 1.0 * Zn_array[indexer(il, in, 2)] )
				/ ( 2.0 * h );

		for (int ix = 1; ix < n_x_pts - 1; ++ix)
		{
			complex<double> fm1 = Zn_array[indexer(il, in, ix-1)];
			complex<double> fp1 = Zn_array[indexer(il, in, ix+1)];
			ddx_Zn_array[indexer(il, in, ix)]
				= ( fp1 - fm1 ) / ( 2.0 * h );
		}

		//approximate first derivative at upper endpoint
		ddx_Zn_array[indexer(il, in, n_x_pts - 1)]
			= ( 3.0 * Zn_array[indexer(il, in, n_x_pts - 1)]
				- 4.0 * Zn_array[indexer(il, in, n_x_pts - 2)]
				+ 1.0 * Zn_array[indexer(il, in, n_x_pts - 3)] )
				/ ( 2.0 * h );

		return;
	}

	void set_d2dx2_Zn(int il, int in)
	{
		double h = x_pts[1] - x_pts[0];

		//approximate first derivative at lower endpoint
		/*d2dx2_Zn_array[indexer(il, in, 0)]
			= ( 1.0 * Zn_array[indexer(il, in, 0)]
				- 2.0 * Zn_array[indexer(il, in, 1)]
				+ 1.0 * Zn_array[indexer(il, in, 2)]
				)
				/ ( h * h );*/
		d2dx2_Zn_array[indexer(il, in, 0)]
			= ( 2.0 * Zn_array[indexer(il, in, 0)]
				- 5.0 * Zn_array[indexer(il, in, 1)]
				+ 4.0 * Zn_array[indexer(il, in, 2)]
				- 1.0 * Zn_array[indexer(il, in, 3)]
				)
				/ ( h * h );

		for (int ix = 1; ix < n_x_pts - 1; ++ix)
		{
			complex<double> fm1 = Zn_array[indexer(il, in, ix-1)];
			complex<double> f0 = Zn_array[indexer(il, in, ix)];
			complex<double> fp1 = Zn_array[indexer(il, in, ix+1)];
			d2dx2_Zn_array[indexer(il, in, ix)]
				= ( fp1 - 2.0 * f0 + fm1 ) / ( h * h );
		}

		//approximate first derivative at upper endpoint
		d2dx2_Zn_array[indexer(il, in, n_x_pts - 1)]
			= ( Zn_array[indexer(il, in, n_x_pts - 1)]
				- 2.0 * Zn_array[indexer(il, in, n_x_pts - 2)]
				+ 1.0 * Zn_array[indexer(il, in, n_x_pts - 3)] )
				/ ( h * h );

		return;
	}

	void set_ddx_x_ddx_Zn(int il, int in)
	{
		for (int ix = 0; ix < n_x_pts; ++ix)
			ddx_x_ddx_Zn_array[indexer(il, in, ix)]
				= x_pts[ix] * d2dx2_Zn_array[indexer(il, in, ix)]
					+ ddx_Zn_array[indexer(il, in, ix)];

		return;
	}


	void set_Zn(int il, int in)
	{
		vector<complex<double> > integrand1(n_x_pts), integrand2(n_x_pts);
		double h = x_pts[1] - x_pts[0];
		for (int ix = 0; ix < n_x_pts; ++ix)
		{
			integrand1[ix] = conj(Zn_array[indexer(il, 0, ix)])
								* ddx_x_ddx_Zn_array[indexer(il, in-1, ix)];
			integrand2[ix] = exp(2.0*i*x_pts[ix])*Zn_array[indexer(il, 0, ix)]
								* ddx_x_ddx_Zn_array[indexer(il, in-1, ix)];
		}

		vector<complex<double> > integrals1(n_x_pts), integrals2(n_x_pts);
		compute_running_integration(&integrand1, &integrals1);
		compute_running_integration(&integrand2, &integrals2);

		for (int ix = 0; ix < n_x_pts; ++ix)
		{
			complex<double> cf1 = -0.5*i*Zn_array[indexer(il, 0, ix)];
			complex<double> cf2 = 0.5*i*exp(-2.0*i*x_pts[ix])*conj(Zn_array[indexer(il, 0, ix)]);

			Zn_array[indexer(il, in, ix)]
				= -cf1*integrals1[ix]-cf2*integrals2[ix];
		}

		return;
	}

	void compute_running_integration(vector<complex<double> > * fpts, vector<complex<double> > * integ_fpts)
	{
		int length = fpts->size();
		double h = x_pts[1] - x_pts[0];
		integ_fpts->resize(length);

		(*integ_fpts)[length - 1] = 0.0;

		for (int ix = length - 2; ix >= 0; --ix)
			(*integ_fpts)[ix]
				= (*integ_fpts)[ix+1] + 0.5*h*((*fpts)[ix] + (*fpts)[ix+1]);

		return;
	}


	void get_all_Zn_and_ddx_Zn(
			vector<complex<double> > * Zn_array_in,
			vector<complex<double> > * ddx_Zn_array_in,
			vector<double> * xpts_in,
			int nmax_in, int lmax_in)
	{
		nmax = nmax_in;
		lmax = lmax_in;
		n_x_pts = xpts_in->size();

		x_pts = *xpts_in;
		Zn_array = vector<complex<double> > ((lmax+1)*nmax*n_x_pts);
		ddx_Zn_array = vector<complex<double> > ((lmax+1)*nmax*n_x_pts);
		d2dx2_Zn_array = vector<complex<double> > ((lmax+1)*nmax*n_x_pts);
		ddx_x_ddx_Zn_array = vector<complex<double> > ((lmax+1)*nmax*n_x_pts);

		//loop over l values eventually
		for (int il = 0; il <= lmax; ++il)
		{
			//set Z0 first
			set_Z0(il);

			//get the rest of the Zn using the recurrence relation
			for (int in = 1; in < nmax; ++in)
			{
				//cout << " * Doing in = " << in << endl;
				set_Zn(il, in);
			}
		}

		/*
		for (int il = 0; il <= lmax; ++il)
		{
			//check against asymptotic expression for F_5(x, x_s)
			double xs = 24.9;
			for (int ix = 0; ix < n_x_pts; ++ix)
			{
				complex<double> sum = 0.0;
		
				for (int in = 0; in < nmax; ++in)
				{
					complex<double> Zn_loc = Zn_array[indexer(il, in, ix)];
					sum += Zn_loc * pow(xs, double(in));
				}
				cout << x_pts[ix] << "   " << sum.real() << "   " << sum.imag() << endl;
			}
		}
		*/

		*Zn_array_in = Zn_array;
		*ddx_Zn_array_in = ddx_Zn_array;

		return;
	}

}

//End of file
