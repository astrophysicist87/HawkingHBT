#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <fenv.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_sf_bessel.h>

#include "chebyshev_library.h"
#include "gauss_quadrature.h"
#include "Zn.h"
#include "Yrt.h"
#include "Xn.h"
#include "main.h"

using namespace std;
using namespace csf;

const complex<double> i(0.0,1.0);
const double xmin = 1.001;
const double xmax = 10.0;
const int nmax = 20;
const int lmax = 0;
const int n_x_pts = 10;
const int rmax = 50;
const int Yrt_nmax = 551;
const int Xn_nmax = 2;
const int smax = Yrt_nmax;
int l = 0;

vector<double> x_pts;
vector<complex<double> > Zn_array((lmax+1)*nmax*n_x_pts);
vector<complex<double> > ddx_Zn_array((lmax+1)*nmax*n_x_pts);
vector<complex<double> > Yrt(rmax*smax*n_x_pts);
vector<complex<double> > Xn(Xn_nmax*n_x_pts);
vector<complex<double> > ddx_Yrt(rmax*smax*n_x_pts);
vector<complex<double> > xinm(rmax*smax);

void set_xinm();

int main (void)
{
	x_pts = vector<double>(n_x_pts);
	Zn_NS::linspace(x_pts, xmin, xmax);

	//Xn_NS::get_Xn(&Xn, &x_pts, Xn_nmax, l);

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

	/*for (int ix = 0; ix < n_x_pts; ++ix)
		cout << x_pts[ix] << "   "
				<< Yrt.at(Yrt_NS::indexer_Yrt(0, 0, ix)).real() << "   "
				<< Yrt.at(Yrt_NS::indexer_Yrt(0, 0, ix)).imag() << endl;*/


	

	double M = 0.5;
	double omega = 1.0;
	double rs = 2.0*M;
	double xs = omega * rs;

	/*
	complex<double> K56 = -2.0*i;
	complex<double> K45 = 0.0;
	for (int n = 0; n < rmax; ++n)
	for (int m = 0; m < smax; ++m)
		K45 += conj(xinm[indexer_xinm(n, m)])
				* pow(xs, double(n))
				* pow(log(xs), double(m));

	double Al = Yrt_NS::Al(l);
	K45 *= pow(-xs, -l) * exp(i * xs) / Al;

	//print out R4
	for (int ix = 0; ix < n_x_pts; ++ix)
	{
		double xpt = x_pts[ix];
		complex<double> sum = 0.0;
		complex<double> prefactor
			= exp(-i * (xpt - xs + xs * log(abs(xpt - xs)))) * pow(-xs, -l) / Al;
		for (int r = 0; r < rmax; ++r)
		for (int t = 0; t < smax; ++t)	//smax is correct
			sum += Yrt.at(Yrt_NS::indexer_Yrt(r, t, ix))
					* pow(xs, r) * pow(log(xs), t);
		sum *= - ( K56 / K45 ) * prefactor;

		cout << xpt << "   " << sum.real() << "   " << sum.imag() << endl;
	}
	*/

	//CHECK: R5 seems to work!!!
	/*
	for (int ix = 0; ix < n_x_pts; ++ix)
	{
		double xpt = x_pts[ix];
		complex<double> sum = 0.0;
		complex<double> prefactor
			= exp(i * (xpt + xs * log(abs(xpt - xs))));

		for (int n = 0; n < nmax; ++n)
			sum += Zn_array.at(Zn_NS::indexer(l, n, ix)) * pow(xs, n);

		sum *= prefactor;

		cout << xpt << "   " << sum.real() << "   " << sum.imag() << endl;
	}
	*/

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
				int nxinf = n_x_pts - 1;
				double x_inf = x_pts[nxinf];
				//cout << n << "   " << m << "   " << n_x_pts - 1 << "   "
				//		<< Yrt_NS::indexer_Yrt(n,m, n_x_pts - 1) << "   " << Yrt.size() << endl;
				xinm[indexer_xinm(n, m)]
					= x_inf*x_inf*
						( Yrt.at(Yrt_NS::indexer_Yrt(n, m, nxinf))
							* conj(ddx_Zn_array.at(Zn_NS::indexer(l, 0, nxinf)))
						- ddx_Yrt.at(Yrt_NS::indexer_Yrt(n, m, nxinf))
							* conj(Zn_array.at(Zn_NS::indexer(l, 0, nxinf))) 
						- 2.0 * i * Yrt.at(Yrt_NS::indexer_Yrt(n, m, nxinf))
							* conj(Zn_array.at(Zn_NS::indexer(l, 0, nxinf)))
						);
				//cout << n << "   " << m << "   "
						/*<< Yrt.at(Yrt_NS::indexer_Yrt(n, m, n_x_pts - 1)) << "   "
						<< conj(ddx_Zn_array.at(Zn_NS::indexer(l, 0, n_x_pts-1))) << "   "
						<< - ddx_Yrt.at(Yrt_NS::indexer_Yrt(n, m, n_x_pts - 1)) << "   "
						<< conj(Zn_array.at(Zn_NS::indexer(l, 0, n_x_pts-1))) << "   "
						<< - 2.0 * i * Yrt.at(Yrt_NS::indexer_Yrt(n, m, n_x_pts-1)) << "   "
						<< conj(Zn_array.at(Zn_NS::indexer(l, 0, n_x_pts-1))) << "   "*/
				//		<< xinm[indexer_xinm(n, m)] << endl;
			//}
	
		}
	}
	return;
}


//End of file
