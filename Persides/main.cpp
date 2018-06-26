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
#include "main.h"

using namespace std;
using namespace csf;

const complex<double> i(0.0,1.0);
const double xmin = 1.0;
const double xmax = 25.0;
const int nmax = 20;
const int lmax = 0;
const int n_x_pts = 10000;
const int rmax = 50;
const int Yrt_nmax = 3;
const int smax = Yrt_nmax;
int l = 0;

vector<double> x_pts;
vector<complex<double> > Zn_array((lmax+1)*nmax*n_x_pts);
vector<complex<double> > ddx_Zn_array((lmax+1)*nmax*n_x_pts);
vector<complex<double> > Yrt(rmax*smax*n_x_pts);
vector<complex<double> > ddx_Yrt(rmax*smax*n_x_pts);
vector<complex<double> > xinm(rmax*smax);

void set_xinm();

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
	
	Yrt_NS::get_Yrt(&Yrt, &ddx_Yrt, &x_pts, Yrt_nmax, rmax, l);

	set_xinm();

	//cout << "Exiting normally." << endl;
	
	return 0;
}

void set_xinm()
{
	//xinm = vector<complex<double> >(rmax*smax);
	for (int n = 0; n < rmax; ++n)
	for (int m = 0; m < smax; ++m)
	{
		if (n<m)	//Eq.(53) of Persides' Paper II
			continue;
		else
		{
			//for (int ix = 0; ix < n_x_pts; ++ix)
			//{
				/*complex<double> sum = 0.0;
				for (int s = 0; s < smax; ++s)
				{
		
				}
				xinm[indexer_xinm(t, m, ix)]
					= conj(sum);*/

				//try simple version: Eq.(54) of Persides' Paper II
				double x_inf = x_pts[n_x_pts-1];
				//cout << n << "   " << m << "   " << n_x_pts - 1 << "   "
				//		<< Yrt_NS::indexer_Yrt(n,m, n_x_pts - 1) << "   " << Yrt.size() << endl;
				xinm[indexer_xinm(n, m)]
					= x_inf*x_inf*
						( Yrt.at(Yrt_NS::indexer_Yrt(n, m, n_x_pts - 1))
							* conj(ddx_Zn_array.at(Zn_NS::indexer(l, 0, n_x_pts-1)))
						- ddx_Yrt.at(Yrt_NS::indexer_Yrt(n, m, n_x_pts - 1))
							* conj(Zn_array.at(Zn_NS::indexer(l, 0, n_x_pts-1))) 
						- 2.0 * i * Yrt.at(Yrt_NS::indexer_Yrt(n, m, n_x_pts-1))
							* conj(Zn_array.at(Zn_NS::indexer(l, 0, n_x_pts-1)))
						);
				cout << n << "   " << m << "   " << xinm[indexer_xinm(n, m)] << endl;
			//}
	
		}
	}
	return;
}


//End of file
