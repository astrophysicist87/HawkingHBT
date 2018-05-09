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
#include "Yrt.h"

using namespace std;
using namespace csf;

const double xmin = 1.0;
const double xmax = 25.0;
const int nmax = 20;
const int lmax = 0;
const int n_x_pts = 10000;
const int rmax = 50;
const int Yrt_nmax = 3;
const int smax = Yrt_nmax;
int l = 2;

vector<double> x_pts;
vector<complex<double> > Zn_array((lmax+1)*nmax*n_x_pts);
vector<complex<double> > ddx_Zn_array((lmax+1)*nmax*n_x_pts);
vector<complex<double> > Yrt(rmax*smax*n_x_pts);
vector<complex<double> > ddx_Yrt(rmax*smax*n_x_pts);

int main (void)
{
	x_pts = vector<double>(n_x_pts);
	Zn_NS::linspace(x_pts, xmin, xmax);

	Zn_NS::get_all_Zn_and_ddx_Zn(&Zn_array, &ddx_Zn_array, &x_pts, nmax, lmax);

	/*
	//loop over l values eventually
	for (int il = 0; il <= lmax; ++il)
	{
		//check against asymptotic expression for F_5(x, x_s)
		double xs = 24.9;
		for (int ix = 0; ix < n_x_pts; ++ix)
		{
			complex<double> sum = 0.0;
		
			for (int in = 0; in < nmax; ++in)
			{
				complex<double> Zn_loc = Zn_array[Zn_NS::indexer(il, in, ix)];
				sum += Zn_loc * pow(xs, double(in));
			}
			cout << x_pts[ix] << "   " << sum.real() << "   " << sum.imag() << endl;
		}
	}
	*/
	
	Yrt_NS::get_Yrt_and_xinm(&Yrt, &ddx_Yrt, &x_pts, Yrt_nmax, rmax, l);

	//cout << "Exiting normally." << endl;
	
	return 0;
}

//End of file
