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
vector<complex<double> > ddx_Zn_array(lmax*nmax*n_x_pts);
vector<complex<double> > ddx_x_ddx_Zn_array(lmax*nmax*n_x_pts);

int main (void)
{
	//set up Chebyshev nodes for evaluating integrals
	//gauss_quadrature(n_x_pts, 5, 0.0, 0.0, 1.0, 20.0*double(n_x_pts)/51.0, x_pts, x_wts);
	gauss_quadrature(n_x_pts, 5, 0.0, 0.0, 1.0, 25.0, x_pts, x_wts);

	//loop over l values eventually
	for (int il = 0; il < lmax; ++il)
	{
		//set Z0 first
		set_Z0(il);

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
					<< setw(16)
					<< setprecision(12)
					<< x_pts[ix] << "   "
					<< Zn_array[indexer(il, 0, ix)].real() << "   "
					<< Zn_array[indexer(il, 0, ix)].imag() << "   "
					//<< ddx_Zn_array[indexer(il, 0, ix)].real() << "   "
					//<< ddx_Zn_array[indexer(il, 0, ix)].imag() << "   "
					//<< ddx_x_ddx_Zn_array[indexer(il, 0, ix)].real() << "   "
					//<< ddx_x_ddx_Zn_array[indexer(il, 0, ix)].imag() 
					<< Z0_deriv1_array[indexer(il, ix)].real() << "   "
					<< Z0_deriv1_array[indexer(il, ix)].imag() << "   "
					<< Z0_deriv2_array[indexer(il, ix)].real() << "   "
					<< Z0_deriv2_array[indexer(il, ix)].imag() << endl;
		*/
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
		Z0_deriv1_array[indexer(l, ix)] = eix * ( coeff1 * h2l + coeff2 * h2lp1 ) / x;
		Z0_deriv2_array[indexer(l, ix)] = eix * ( coeff1 * h1l + coeff2 * h1lp1 ) / x;
	}

	//set first derivative
	//set_ddx_Zn(l, 0);

	//set second derivative combination
	//set_ddx_x_ddx_Zn(l, 0);

	return;
}

/*void set_ddx_Zn(int il, int in)
{
	ddx_Zn_array[indexer(il, in, 0)]
		= ( Zn_array[indexer(il, in, 1)]
			- Zn_array[indexer(il, in, 0)] )
			/ ( x_pts[1] - x_pts[0] );
	for (int ix = 1; ix < n_x_pts; ++ix)
			ddx_Zn_array[indexer(il, in, ix)]
				= ( Zn_array[indexer(il, in, ix)]
					- Zn_array[indexer(il, in, ix-1)] )
					/ ( x_pts[ix] - x_pts[ix-1] );

	return;
}

void set_ddx_x_ddx_Zn(int il, int in)
{
	ddx_x_ddx_Zn_array[indexer(il, in, 0)]
		= ( x_pts[1]*ddx_Zn_array[indexer(il, in, 1)]
			- x_pts[0]*ddx_Zn_array[indexer(il, in, 0)] )
			/ ( x_pts[1] - x_pts[0] );
	for (int ix = 1; ix < n_x_pts; ++ix)
			ddx_x_ddx_Zn_array[indexer(il, in, ix)]
				= ( x_pts[ix]*ddx_Zn_array[indexer(il, in, ix)]
					- x_pts[ix-1]*ddx_Zn_array[indexer(il, in, ix-1)] )
					/ ( x_pts[ix] - x_pts[ix-1] );

	return;
}*/

void set_Zn(int il, int in)
{
	for (int ix = 0; ix < n_x_pts; ++ix)
	{
		double x_loc = x_pts[ix];
		complex<double> Z0_loc = Zn_array[indexer(il, 0, ix)];
		complex<double> term1 = x_loc * Z0_loc * conj(Z0_loc) * Zn_array[indexer(il, in - 1, ix)];
		complex<double> term2 = 0.0, term3 = 0.0;
		for (int ixp = ix+1; ixp < n_x_pts; ++ixp)
		{
			term2 += x_wts[ixp] * Z0_deriv1_array[indexer(il, ixp)] * Zn_array[indexer(il, in - 1, ixp)];
			term3 += x_wts[ixp] * Z0_deriv2_array[indexer(il, ixp)] * Zn_array[indexer(il, in - 1, ixp)];
		}
		term2 *= -0.5*i*Z0_loc;
		term3 *= 0.5*i*exp(-2.0*i*x_loc)*conj(Z0_loc);

		//note minus signs: integration limits reversed!
		Zn_array[indexer(il, in, ix)] = term1 - term2 - term3;
	}

	return;
}


