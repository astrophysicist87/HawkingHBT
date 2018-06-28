#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sum.h>

#include "Xn.h"

using namespace std;

#define USE_SERIES_ACCELERATION		true

namespace Xn_NS
{
	int nmax, l, n_x_pts;
	double sigma_l = 0.0;

	vector<double> x_pts;

	vector<complex<double> > Xn;
	vector<double> Xtn;		//hopefully keeps memory management a bit simpler

	vector<double> Pl_arg, Ql_arg;

	void set_Pl_arg_and_Ql_arg(int l_local)
	{
		//set coefficients for Pl_arg evaluation
		double ppr[l_local+1];
		for (int il = 0; il <= l_local; ++il)
			ppr[il] = pp(il, l_local);

		//now set Pl_arg
		for (int ix = 0; ix < n_x_pts; ++ix)
			Pl_arg[ix] = eval_poly(ppr, l_local, x_pts[ix]);

		//now set Ql_arg
		double prefactor = pow(-1.0, l_local+1);
		for (int ix = 0; ix < n_x_pts; ++ix)
			Ql_arg[ix] = prefactor * gsl_sf_legendre_Ql(l_local, 2.0*x_pts[ix]-1.0);

		/*
		//check
		for (int ix = 0; ix < n_x_pts; ++ix)
			cout << l_local << "   " << x_pts[ix] << "   " << Pl_arg[ix] << "   " << Ql_arg[ix] << endl;
		//if (1) exit (8);
		*/

		return;
	}

	void set_Xt0_and_Xt1(double sigma_l)
	{
		//set Xt0
		double prefactor = pow(-1.0, l);
		for (int ix = 0; ix < n_x_pts; ++ix)
			Xtn[indexer_Xtn(0, ix)] = prefactor * Pl_arg[ix];

		//set Xt1
		for (int ix = 0; ix < n_x_pts; ++ix)
		{
			double xpt = x_pts[ix];
			Xtn[indexer_Xtn(1, ix)]
				= prefactor
					* ( ( 1.0 - 2.0*sigma_l - xpt - log(xpt - 1.0) )
						* Pl_arg[ix] + 2.0 * Ql_arg[ix] );
		}

		for (int ix = 0; ix < n_x_pts; ++ix)
			cout << x_pts[ix] << "   "
					<< Xtn[indexer_Xtn(0, ix)] << "   "
					<< Xtn[indexer_Xtn(1, ix)] << endl;
		if (1) exit (8);

		return;
	}

	/*void set_Xtn(int n)
	{
		for (int ix = 0; ix < n_x_pts; ++ix)
		{

		}

		return;
	}

	//use Xtn and rescale to get Xn
	void set_Xn(int n)
	{
		for (int ix = 0; ix < n_x_pts; ++ix)
		{

		}

		return;
	}*/

	void get_Xn(
			vector<complex<double> > * Xn_in,
			vector<double> * x_pts_in,
			int nmax_in, int l_in)
	{
		nmax = nmax_in;
		l = l_in;
		n_x_pts = x_pts_in->size();

		x_pts = *x_pts_in; 

		Xn = vector<complex<double> >(nmax*n_x_pts);
		Xtn = vector<double>(nmax*n_x_pts);

		Pl_arg = vector<double>(n_x_pts);
		Ql_arg = vector<double>(n_x_pts);

		//testing
		//int lmax = 10;
		//for (int chosen_l = 0; chosen_l <= lmax; ++chosen_l)
		//	set_Pl_arg_and_Ql_arg(chosen_l);
		set_Pl_arg_and_Ql_arg(l);

		for (int m = 1; m <= l; ++m)
			sigma_l += 1.0/m;
		//cout << sigma_l << endl;

		//cout << "Start." << endl;
		set_Xt0_and_Xt1(sigma_l);
	
		/*
		//cout << "Next." << endl;
		for (int n = 2; n < nmax; ++n)
		{
			cout << "n = " << n << endl;
			set_Xn(n);
		}
		//cout << "Here." << endl;

		//finally, be sure to copy over results!
		for (int iXn = 0; iXn < Xn.size(); ++iXn)
			Xn_in->at(iXn) = Xn.at(iXn);
		*/

		return;
	}
}

