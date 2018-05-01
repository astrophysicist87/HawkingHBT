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

vector<double> x_pts, x_wts;
vector<complex<double> > Zn_array(lmax*nmax*n_x_pts);
vector<complex<double> > ddx_Zn_array(lmax*nmax*n_x_pts);
vector<complex<double> > ddx_x_ddx_Zn_array(lmax*nmax*n_x_pts);

int main (void)
{
	gsl_cheb_series * cs;

	//set up Chebyshev nodes for evaluating integrals
	gauss_quadrature(n_x_pts, 1, 0.0, 0.0, xmax, xmin, x_pts, x_wts);
	
	//set up Chebyshev stuff
	//set_up_Chebyshev(n_x_pts);

	//loop over l values eventually
	for (int il = 0; il < lmax; ++il)
	{
		//cout << "Working on il = " << il << endl;
		cs = gsl_cheb_alloc(order);

		//set Z0 first
		set_Z0(il);

		//get the rest of the Zn using the recurrence relation
		for (int in = 1; in < nmax; ++in)
		{
			//cout << " * Doing in = " << in << endl;
			set_Zn(il, in);
		}

		//output/return Zn...
		/*for (int in = 0; in < nmax; ++in)
		for (int ix = 0; ix < n_x_pts; ++ix)
		{
			double re = Zn_array[indexer(il, in, ix)].real();
			double im = Zn_array[indexer(il, in, ix)].imag();
			cout
				<< il << "   "
				<< in << "   "
				<< setw(20)
				<< setprecision(16)
				<< x_pts[ix] << "   "
				<< re << "   "
				<< im << endl;
		}*/

		for (int ix = 0; ix < n_x_pts; ++ix)
			cout << il << "   "
					<< setw(16)
					<< setprecision(12)
					<< x_pts[ix] << "   "
					<< Zn_array[indexer(il, 0, ix)].real() << "   "
					<< Zn_array[indexer(il, 0, ix)].imag() << "   "
					<< ddx_Zn_array[indexer(il, 0, ix)].real() << "   "
					<< ddx_Zn_array[indexer(il, 0, ix)].imag() << "   "
					<< ddx_x_ddx_Zn_array[indexer(il, 0, ix)].real() << "   "
					<< ddx_x_ddx_Zn_array[indexer(il, 0, ix)].imag() << endl;

		gsl_cheb_free(cs);
	}

	//cout << "Exiting normally." << endl;
	
	return 0;
}



//////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////

/*inline int indexer(int il, int ix)
{
	return( il * n_x_pts + ix );
}*/

inline int indexer(int il, int in, int ix)
{
	return( ( il * nmax + in ) * n_x_pts + ix );
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

void convert_array_to_Chebyshev_coeffs(double a_in, double b_in, vector<double> values_on_nodes, gsl_cheb_series * coefficients)
{
	int n_coeffs = (int)values_on_nodes.size();
	//cout << "n_coeffs = " << n_coeffs << endl;
	double coeffs_array[n_coeffs];
	double * local_cheb_coeffs = new double [n_coeffs];

	double default_nodes[n_coeffs];
	for (int k = 0; k < n_coeffs; ++k)
		default_nodes[k] = - cos( M_PI*(2.*(k+1.) - 1.) / (2.*n_coeffs) );

	double nums[n_coeffs*n_coeffs];
	double dens[n_coeffs];

	for (int j = 0; j < n_coeffs; ++j)
	{
		dens[j] = 0.0;
		for (int k = 0; k < n_coeffs; ++k)
		{
			double Tjk = csf::Tfun(j, default_nodes[k]);
			dens[j] += Tjk*Tjk;
			nums[j*n_coeffs+k] = Tjk;
		}
	}

	//////////////////////////////////
	//separate out 0th coefficient for additional factor of 2.0
	coeffs_array[0] = 0.0;
	for (int k = 0; k < n_coeffs; ++k)
		coeffs_array[0] += 2.0 * values_on_nodes[k] * nums[0*n_coeffs+k];

	local_cheb_coeffs[0] = coeffs_array[0] / dens[0];

	for (int j = 1; j < n_coeffs; ++j)
	{
		coeffs_array[j] = 0.0;
		for (int k = 0; k < n_coeffs; ++k)
			coeffs_array[j] += values_on_nodes[k] * nums[j*n_coeffs+k];
		local_cheb_coeffs[j] = coeffs_array[j] / dens[j];
	}

	//set gsl data for chebyshev evaluation
	coefficients->a = a_in;
	coefficients->b = b_in;
	for (int j = 0; j < n_coeffs; ++j)
		coefficients->c[j] = local_cheb_coeffs[j];

	delete [] local_cheb_coeffs;

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
	}

	//set first derivative
	set_ddx_Zn(l, 0);

	//set second derivative combination
	set_ddx_x_ddx_Zn(l, 0);

	return;
}

void set_ddx_Zn(int il, int in)
{
	gsl_cheb_series * cs_Zn_re = gsl_cheb_alloc(order);
	gsl_cheb_series * cs_Zn_im = gsl_cheb_alloc(order);
	gsl_cheb_series * cs_ddx_Zn_re = gsl_cheb_alloc(order);
	gsl_cheb_series * cs_ddx_Zn_im = gsl_cheb_alloc(order);

	vector<double> Zn_realPart(n_x_pts), Zn_imagPart(n_x_pts);

	for (int ix = 0; ix < n_x_pts; ++ix)
	{
		complex<double> Zn_loc = Zn_array[indexer(il, in, ix)];
		Zn_realPart[ix] = Zn_loc.real();
		Zn_imagPart[ix] = Zn_loc.imag();
	}
	
	//get corresponding chebyshev coefficients
	convert_array_to_Chebyshev_coeffs(xmax, xmin, Zn_realPart, cs_Zn_re);
	convert_array_to_Chebyshev_coeffs(xmax, xmin, Zn_imagPart, cs_Zn_im);

	//take derivative
	gsl_cheb_calc_deriv(cs_ddx_Zn_re, cs_Zn_re);
	gsl_cheb_calc_deriv(cs_ddx_Zn_im, cs_Zn_im);

	for (int ix = 0; ix < n_x_pts; ++ix)
	{
		double x_loc = x_pts[ix];
		double ddx_Zn_re_loc = gsl_cheb_eval(cs_ddx_Zn_re, x_loc);
		double ddx_Zn_im_loc = gsl_cheb_eval(cs_ddx_Zn_im, x_loc);
		ddx_Zn_array[indexer(il, in, ix)]
			= complex<double>(ddx_Zn_re_loc, ddx_Zn_im_loc);
	}

	gsl_cheb_free(cs_Zn_re);
	gsl_cheb_free(cs_Zn_im);
	gsl_cheb_free(cs_ddx_Zn_re);
	gsl_cheb_free(cs_ddx_Zn_im);

	return;
}

void set_ddx_x_ddx_Zn(int il, int in)
{
	gsl_cheb_series * cs_x_ddx_Zn_re = gsl_cheb_alloc(order);
	gsl_cheb_series * cs_x_ddx_Zn_im = gsl_cheb_alloc(order);
	gsl_cheb_series * cs_ddx_x_ddx_Zn_re = gsl_cheb_alloc(order);
	gsl_cheb_series * cs_ddx_x_ddx_Zn_im = gsl_cheb_alloc(order);

	vector<double> x_ddx_Zn_realPart(n_x_pts), x_ddx_Zn_imagPart(n_x_pts);

	for (int ix = 0; ix < n_x_pts; ++ix)
	{
		double x_loc = x_pts[ix];
		complex<double> x_ddx_Zn_loc = x_loc*ddx_Zn_array[indexer(il, in, ix)];
		x_ddx_Zn_realPart[ix] = x_ddx_Zn_loc.real();
		x_ddx_Zn_imagPart[ix] = x_ddx_Zn_loc.imag();
	}
	
	//get corresponding chebyshev coefficients
	convert_array_to_Chebyshev_coeffs(xmax, xmin, x_ddx_Zn_realPart, cs_x_ddx_Zn_re);
	convert_array_to_Chebyshev_coeffs(xmax, xmin, x_ddx_Zn_imagPart, cs_x_ddx_Zn_im);

	//take derivative
	gsl_cheb_calc_deriv(cs_ddx_x_ddx_Zn_re, cs_x_ddx_Zn_re);
	gsl_cheb_calc_deriv(cs_ddx_x_ddx_Zn_im, cs_x_ddx_Zn_im);

	for (int ix = 0; ix < n_x_pts; ++ix)
	{
		double x_loc = x_pts[ix];
		double ddx_x_ddx_Zn_re_loc = gsl_cheb_eval(cs_ddx_x_ddx_Zn_re, x_loc);
		double ddx_x_ddx_Zn_im_loc = gsl_cheb_eval(cs_ddx_x_ddx_Zn_im, x_loc);
		ddx_x_ddx_Zn_array[indexer(il, in, ix)]
			= complex<double>(ddx_x_ddx_Zn_re_loc, ddx_x_ddx_Zn_im_loc);
	}

	gsl_cheb_free(cs_x_ddx_Zn_re);
	gsl_cheb_free(cs_x_ddx_Zn_im);
	gsl_cheb_free(cs_ddx_x_ddx_Zn_re);
	gsl_cheb_free(cs_ddx_x_ddx_Zn_im);

	return;
}


void set_integrand_1_coeffs(int il, int in,
							gsl_cheb_series * csIntegrand_re,
							gsl_cheb_series * csIntegrand_im)
{
	//set this integrand, real and imaginary parts
	vector<double> integrand_realPart(n_x_pts), integrand_imagPart(n_x_pts);
	for (int ix = 0; ix < n_x_pts; ++ix)
	{
		complex<double> integrand
			= conj(Zn_array[indexer(il, 0, ix)])
				* ddx_x_ddx_Zn_array[indexer(il, in-1, ix)];
		integrand_realPart[ix] = integrand.real();
		integrand_imagPart[ix] = integrand.imag();
	}

	//convert real and imaginary parts to corresponding Chebyshev coefficient arrays
	convert_array_to_Chebyshev_coeffs(xmax, xmin, integrand_realPart, csIntegrand_re);
	convert_array_to_Chebyshev_coeffs(xmax, xmin, integrand_imagPart, csIntegrand_im);

	return;
}

void set_integrand_2_coeffs(int il, int in,
							gsl_cheb_series * csIntegrand_re,
							gsl_cheb_series * csIntegrand_im)
{
	//set this integrand, real and imaginary parts
	vector<double> integrand_realPart(n_x_pts), integrand_imagPart(n_x_pts);
	for (int ix = 0; ix < n_x_pts; ++ix)
	{
		double x_loc = x_pts[ix];
		complex<double> integrand
			= exp(2.0*i*x_loc) * Zn_array[indexer(il, 0, ix)]
				* ddx_x_ddx_Zn_array[indexer(il, in-1, ix)];
		integrand_realPart[ix] = integrand.real();
		integrand_imagPart[ix] = integrand.imag();
	}

	//convert real and imaginary parts to corresponding Chebyshev coefficient arrays
	convert_array_to_Chebyshev_coeffs(xmax, xmin, integrand_realPart, csIntegrand_re);
	convert_array_to_Chebyshev_coeffs(xmax, xmin, integrand_imagPart, csIntegrand_im);

	return;
}

void set_Zn(int l, int in)
{
	//initialize chebyshev coefficient arrays
	gsl_cheb_series * csIntegrand1_re = gsl_cheb_alloc(order);
	gsl_cheb_series * csIntegrand1_im = gsl_cheb_alloc(order);
	gsl_cheb_series * csIntegrand2_re = gsl_cheb_alloc(order);
	gsl_cheb_series * csIntegrand2_im = gsl_cheb_alloc(order);
	gsl_cheb_series * csIntegral1_re = gsl_cheb_alloc(order);
	gsl_cheb_series * csIntegral1_im = gsl_cheb_alloc(order);
	gsl_cheb_series * csIntegral2_re = gsl_cheb_alloc(order);
	gsl_cheb_series * csIntegral2_im = gsl_cheb_alloc(order);

	//set integrand #1 (real and imaginary parts)
	set_integrand_1_coeffs(l, in,
							csIntegrand1_re,
							csIntegrand1_im);

	//set integrand #2 (real and imaginary parts)
	set_integrand_2_coeffs(l, in,
							csIntegrand2_re,
							csIntegrand2_im);

	//get integral coefficients
	gsl_cheb_calc_integ(csIntegral1_re, csIntegrand1_re);
	gsl_cheb_calc_integ(csIntegral1_im, csIntegrand1_im);
	gsl_cheb_calc_integ(csIntegral2_re, csIntegrand2_re);
	gsl_cheb_calc_integ(csIntegral2_im, csIntegrand2_im);

	//do it for each point x
	for (int ix = 0; ix < n_x_pts; ++ix)
	{
		double x = x_pts[ix];

		//compute integral #1
		double integral1_re = gsl_cheb_eval(csIntegral1_re, x);
		double integral1_im = gsl_cheb_eval(csIntegral1_im, x);
		complex<double> integral1(integral1_re, integral1_im);

		//compute integral #2
		double integral2_re = gsl_cheb_eval(csIntegral2_re, x);
		double integral2_im = gsl_cheb_eval(csIntegral2_im, x);
		complex<double> integral2(integral2_re, integral2_im);

		//assign result to appropriate Zn_array cell
		complex<double> Z0 = Zn_array[indexer(l, 0, ix)];
		complex<double> Z0c = conj(Z0);
		Zn_array[indexer(l, in, ix)] =
				- 0.5 * i * Z0 * integral1
				+ 0.5 * i * Z0c * exp(-2.0*i*x) * integral2;
	}

	//set first derivative
	set_ddx_Zn(l, in);

	//set second derivative combination
	set_ddx_x_ddx_Zn(l, in);

	gsl_cheb_free(csIntegrand1_re);
	gsl_cheb_free(csIntegrand1_im);
	gsl_cheb_free(csIntegrand2_re);
	gsl_cheb_free(csIntegrand2_im);
	gsl_cheb_free(csIntegral1_re);
	gsl_cheb_free(csIntegral1_im);
	gsl_cheb_free(csIntegral2_re);
	gsl_cheb_free(csIntegral2_im);

	return;
}



