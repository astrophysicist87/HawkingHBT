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
vector<complex<double> > Z0_deriv1_array(lmax*n_x_pts);
vector<complex<double> > Z0_deriv2_array(lmax*n_x_pts);

int main (void)
{
	gsl_cheb_series * cs;

	//set up Chebyshev nodes for evaluating integrals
	gauss_quadrature(n_x_pts, 1, 0.0, 0.0, xmax, 1.0, x_pts, x_wts);
	
	//set up Chebyshev stuff
	//set_up_Chebyshev(n_x_pts);

	//loop over l values eventually
	for (int il = 0; il < lmax; ++il)
	{
		//cout << "Working on il = " << il << endl;
		cs = gsl_cheb_alloc(order);

		//set Z0 first
		set_Z0_and_derivatives(il);

		//get the rest of the Zn using the recurrence relation
		for (int in = 1; in < nmax; ++in)
		{
			//cout << " * Doing in = " << in << endl;
			set_Zn(il, in);
		}

		//output/return Zn...
		for (int in = 0; in < nmax; ++in)
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
		}

		/*for (int ix = 0; ix < n_x_pts; ++ix)
			cout << il << "   "
					<< x_pts[ix] << "   "
					<< Z0_deriv1_array[indexer(il, ix)].real() << "   "
					<< Z0_deriv1_array[indexer(il, ix)].imag() << "   "
					<< Z0_deriv2_array[indexer(il, ix)].real() << "   "
					<< Z0_deriv2_array[indexer(il, ix)].imag() << endl;*/

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

inline int indexer(int il, int ix)
{
	return( il * n_x_pts + ix );
}

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

void convert_array_to_Chebyshev_coeffs(vector<double> values_on_nodes, gsl_cheb_series * coefficients)
{
	int n_coeffs = (int)values_on_nodes.size();
	//cout << "n_coeffs = " << n_coeffs << endl;
	double coeffs_array[n_coeffs];
	double * local_cheb_coeffs = new double [n_coeffs];

	double default_nodes[n_coeffs];
	for (int k = 0; k < n_coeffs; ++k)
	{
		default_nodes[k] = - cos( M_PI*(2.*(k+1.) - 1.) / (2.*n_coeffs) );
		//these would be the node values, but they're not needed here...
		//alpha_pts[k] = 0.5*(default_nodes[k] + 1.0)*(alpha_max - alpha_min) + alpha_min;
	}
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
	//coefficients->a = 1.0;
	//coefficients->b = xmax;
	coefficients->a = xmax;
	coefficients->b = 1.0;
	//coefficients->c = local_cheb_coeffs;
	for (int j = 0; j < n_coeffs; ++j)
		coefficients->c[j] = local_cheb_coeffs[j];

	delete [] local_cheb_coeffs;

	return;
}

void set_Z0_and_derivatives(int l)
{
	//cout << "\t - Made it to set_Z0_and_derivatives()" << endl;
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
		Z0_deriv1_array[indexer(l, ix)] = eix * ( coeff1 * h2l + coeff2 * h2lp1 ) / x;
		Z0_deriv2_array[indexer(l, ix)] = eix * ( coeff1 * h1l + coeff2 * h1lp1 ) / x;
	}
	//cout << "\t - Made out of set_Z0_and_derivatives()" << endl;

	return;
}

void set_integrand_coeffs(vector<complex<double> > Z0_deriv_array,
							int il, int in,
							gsl_cheb_series * csIntegrand_re,
							gsl_cheb_series * csIntegrand_im)
{
	//set this integrand, real and imaginary parts
	vector<double> integrand_realPart(n_x_pts), integrand_imagPart(n_x_pts);
	for (int ix = 0; ix < n_x_pts; ++ix)
	{
		complex<double> integrand
			= Z0_deriv_array[indexer(il, ix)]
				* Zn_array[indexer(il, in-1, ix)];
		integrand_realPart[ix] = integrand.real();
		integrand_imagPart[ix] = integrand.imag();
//cout << "integrand points: " << il << "   " << in << "   " << x_pts[ix] << "   "
//		<< integrand.real() << "   " << integrand.imag() << endl;
	}

	//convert real and imaginary parts to corresponding Chebyshev coefficient arrays
	convert_array_to_Chebyshev_coeffs(integrand_realPart, csIntegrand_re);
	convert_array_to_Chebyshev_coeffs(integrand_imagPart, csIntegrand_im);

	//check coeffs calculation
	for (int ix = 0; ix < 10000; ++ix)
	{
		double hw = 0.5*(1.0-xmax);
		double cen = 0.5*(1.0+xmax);
//cout << "integrand interp points: " << il << "   " << in << "   " << cen+0.0001*hw*ix << "   "
//		<< gsl_cheb_eval(csIntegrand_re, cen+0.0001*hw*ix) << "   "
//		<< gsl_cheb_eval(csIntegrand_im, cen+0.0001*hw*ix) << endl;
		
	}


	return;
}

void set_Zn(int l, int in)
{
	//cout << "\t - Made it to set_Zn()" << endl;
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
	set_integrand_coeffs(Z0_deriv1_array, l, in,
							csIntegrand1_re,
							csIntegrand1_im);

	//set integrand #2 (real and imaginary parts)
	set_integrand_coeffs(Z0_deriv2_array, l, in,
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
		Zn_array[indexer(l, in, ix)]
			= x * Z0 * Z0c * Zn_array[indexer(l, in-1, ix)]
				- 0.5 * i * Z0 * integral1
				+ 0.5 * i * Z0c * exp(-2.0*i*x) * integral2;
	}

	gsl_cheb_free(csIntegrand1_re);
	gsl_cheb_free(csIntegrand1_im);
	gsl_cheb_free(csIntegrand2_re);
	gsl_cheb_free(csIntegrand2_im);
	gsl_cheb_free(csIntegral1_re);
	gsl_cheb_free(csIntegral1_im);
	gsl_cheb_free(csIntegral2_re);
	gsl_cheb_free(csIntegral2_im);

	//cout << "\t - Made out of set_Zn()" << endl;

	return;
}



