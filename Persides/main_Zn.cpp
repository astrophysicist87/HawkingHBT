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
//using namespace Zn_NS;

vector<double> x_pts, x_wts;
vector<complex<double> > Zn_array((lmax+1)*nmax*n_x_pts);
vector<complex<double> > Z0_deriv1_array((lmax+1)*n_x_pts);
vector<complex<double> > Z0_deriv2_array((lmax+1)*n_x_pts);
vector<complex<double> > ddx_Zn_array((lmax+1)*nmax*n_x_pts);
vector<complex<double> > d2dx2_Zn_array((lmax+1)*nmax*n_x_pts);
vector<complex<double> > ddx_x_ddx_Zn_array((lmax+1)*nmax*n_x_pts);

int main (void)
{
	x_pts = vector<double>(n_x_pts);
	linspace(x_pts, xmin, xmax);



	
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

	}


	//loop over l values eventually
	for (int il = 0; il <= lmax; ++il)
	{

		//check against asymptotic expression for F_5(x, x_s)
		///*
		double xs = 24.9;
		for (int ix = 0; ix < n_x_pts; ++ix)
		{
			complex<double> sum = 0.0;
		
			for (int in = 0; in < nmax; ++in)
			{
				complex<double> Zn_loc = Zn_NS::Zn_array[indexer(il, in, ix)];
				sum += Zn_loc * pow(xs, double(in));
			}
			cout << x_pts[ix] << "   " << sum.real() << "   " << sum.imag() << endl;
		}
		//*/

	}

	//cout << "Exiting normally." << endl;
	
	return 0;
}

//End of file
