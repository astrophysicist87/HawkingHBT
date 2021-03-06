#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <fenv.h>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sum.h>

#include "Yrt.h"

using namespace std;

#define USE_SERIES_ACCELERATION		true

namespace Yrt_NS
{
	int rmax, nmax, smax, tmax, l, n_x_pts;
	double sigma_l = 0.0;

	vector<double> x_pts;

	vector<complex<double> > B1D, C1D;
	vector<complex<double> > A3D, B3D, C3D;
	//vector<complex<double> > A4D;

	vector<complex<double> > Yrt, ddx_Yrt;

	void initialize_A3D_starting_values()
	{
		//set n==0,s==0
		for (int r = 0; r < rmax; ++r)
			A3D.at(indexer_A3D(0, r, 0))
				= pow(-1.0, l)*pp(r);

		//set n==1,s==0
		for (int r = 0; r < rmax; ++r)
		{
			double small_sum = 0.0;
			for (int m = 0; m <= r; ++m)
				small_sum += pp(r-m-2)/(m+1.0);
		
			A3D.at(indexer_A3D(1, r, 0))
				= (1.0-2.0*sigma_l)*i*pow(-1.0,l)*pp(r-1)
					+ 2.0*i*pow(-1.0,l)*qp(r-2*l-2)
					- i*pow(-1.0,l)*pp(r)
					+ i*pow(-1.0,l)*small_sum;
		}

		//set n==1,s==1
		for (int r = 0; r < rmax; ++r)
			A3D.at(indexer_A3D(1, r, 1))
				= i*pow(-1.0, l+1)*pp(r-1);

		/*
		for (int n = 0; n < nmax; ++n)
		for (int r = 0; r < rmax; ++r)
		for (int s = 0; s < smax; ++s)
			cout << "A3D[" << n << "," << r << "," << s << "] = " << A3D.at(indexer_A3D(n, r, s)) << endl;
		*/

		return;
	}

	void set_B1D(int n)
	{
		complex<double> sum = 0.0;
		for (int u = 0; u < rmax; ++u)
		{
			int vupper = min(min(n-1, u), smax-1);
			for (int t = 0; t <= u+l; ++t)	//double check
			//for (int v = 0; v < smax; ++v)
			for (int v = 0; v <= vupper; ++v)
			{
				//cout << "B1D: " << (l-t+u+1.0) << "   " << pp(t-u) << "   "
				//		<< A3D.at(indexer_A3D(n-1,u,v)) << "   " << D3D(2*l+n-t,v,0) << endl;

				sum += (l-t+u+1.0)*pp(t-u)*A3D.at(indexer_A3D(n-1,u,v))*D3D(2*l+n-t,v,0);
			}
		}

		B1D.at(n) = sum;
		//cout << "B1D[" << n << "] = " << B1D.at(n) << endl;

		return;
	}

	void compute_levin_sum(vector<double> computed_terms, double * sum_accel, double * err)
	{
		//cout << "In compute_levin_sum(...): N = "<< computed_terms.size() << endl;
		const int N = computed_terms.size();
		double t[N];
		double sum = 0;

		gsl_sum_levin_u_workspace * w = gsl_sum_levin_u_alloc (N);

		for (int n = 0; n < N; n++)
		{
			t[n] = computed_terms[n];
			sum += t[n];
		}

		gsl_sum_levin_u_accel (t, N, w, sum_accel, err);

		//printf ("term-by-term sum = % .16f using %d terms\n", sum, N);

		//printf ("term-by-term sum = % .16f using %zu terms\n", w->sum_plain, w->terms_used);

		//printf ("accelerated sum  = % .16f using %zu terms\n", *sum_accel, w->terms_used);

		//printf ("estimated error  = % .16f\n", *err);

		gsl_sum_levin_u_free (w);

		//cout << "Exiting compute_levin_sum(...)." << endl;

		return;
	}

	void set_C1D(int n)
	{
		complex<double> sum = 0.0, preliminary_sum = 0.0, sum_accel = 0.0;
		vector<double> running_t_terms;
		for (int t = 0; t < tmax; ++t) //double check
		{
			complex<double> sum_t = 0.0;
			//for (int u = 0; u < rmax; ++u)
			for (int u = 0; u <= t; ++u)
			{
				//complex<double> sum_ut = 0.0;
				int vupper = min(min(n-1, u), smax-1);
				//for (int v = 0; v < smax; ++v)
				for (int v = 0; v <= vupper; ++v)
				{
					//cout << "C1D: " << (-l-t+u) << "   " << qp(t-u) << "   "
					//		<< A3D.at(indexer_A3D(n-1,u,v)) << "   " << D3D(n-1-t,v,0) << endl;
					//sum += (-l-t+u)*qp(t-u)*A3D.at(indexer_A3D(n-1,u,v))*D3D(n-1-t,v,0);
					sum_t += (-l-t+u)*qp(t-u)*A3D.at(indexer_A3D(n-1,u,v))*D3D(n-1-t,v,0);
					//sum_ut += (-l-t+u)*qp(t-u)*A3D.at(indexer_A3D(n-1,u,v))*D3D(n-1-t,v,0);
				}

			}
			if ( USE_SERIES_ACCELERATION and t > n )
			{
				running_t_terms.push_back( (pow(i,n+1.0)*sum_t).real() );
			}
			else
				preliminary_sum += sum_t;
		}

		double sum_accel_im = 0.0, err = 0.0;
		if ( USE_SERIES_ACCELERATION and running_t_terms.size() > 0 )
			compute_levin_sum(running_t_terms, &sum_accel_im, &err);

		C1D.at(n) = preliminary_sum + pow(i,-n-1.0)*sum_accel_im;
		//cout << "C1D[" << n << "] = "
		//		<< preliminary_sum + pow(i,-n-1.0)*sum_accel_im << endl;

		return;
	}

	void set_A3D(int n, int r, int s)
	{
		complex<double> sum = 0.0;
		for (int t = 0; t <= r; ++t) //double check
			sum += qp(r-t)*B3D.at(indexer_B3D(n, t, s))
					- pp(r-t)*C3D.at(indexer_C3D(n, t, s));

		A3D.at(indexer_A3D(n,r,s)) = 4.0*i*sum;
		cout << "A3D[" << n << "," << r << "," << s << "] = " << A3D.at(indexer_A3D(n,r,s)) << endl;
		return;
	}

	void set_B3D(int n, int t, int s)
	{
		complex<double> sum = 0.0;
		for (int u = 0; u < rmax; ++u)
		{
			int vupper = min(min(n-1, u), smax-1);
			//for (int v = 0; v < smax; ++v)
			for (int v = 0; v <= vupper; ++v)
				sum += (l-t+u+1.0)*pp(t-u)*A3D.at(indexer_A3D(n-1,u,v))*D3D(2*l+n-t,v,s);
		}

		B3D.at(indexer_B3D(n, t, s)) = sum - B1D.at(n)*double(s==0)*double(2*l+n+1==t);
		//cout << "B3D[" << n << "," << t << "," << s << "] = " << B3D.at(indexer_B3D(n, t, s)) << endl;

		return;
	}

	void set_C3D(int n, int t, int s)
	{
		complex<double> sum = 0.0;
		for (int u = 0; u < rmax; ++u)
		{
			int vupper = min(min(n-1, u), smax-1);
			//for (int v = 0; v < smax; ++v)
			for (int v = 0; v <= vupper; ++v)
				sum += (-l-t+u)*qp(t-u)*A3D.at(indexer_A3D(n-1,u,v))*D3D(n-1-t,v,s);
		}

		C3D.at(indexer_C3D(n, t, s)) = sum - C1D.at(n)*double(s==0)*double(n==t);
		//cout << "C3D[" << n << "," << t << "," << s << "] = " << C3D.at(indexer_C3D(n, t, s)) << endl;

		return;
	}


	///*
	void set_Yrt()
	{
		Yrt = vector<complex<double> >(rmax*smax*n_x_pts);
		double prefactor = pow(-1.0,l)*Al(l);
		for (int ix = 0; ix < n_x_pts; ++ix)
		for (int r = 0; r < rmax; ++r)
		for (int t = 0; t < smax; ++t)	//smax is correct
		{
			if (r<t)	//Eq.(47) of Persides' Paper II
				continue;
			double x_loc = x_pts[ix];
			double ln_x_loc = log(x_loc);
			complex<double> sum = 0.0;
			for (int n = 0; n < nmax; ++n)
			for (int s = 0; s < smax; ++s)
			{
				if (t>s)
					continue;
				complex<double> Anrst
					= pow(-1.0,t)
						* gsl_sf_choose(s, t)
						* A3D[indexer_A3D(n,r,s)];
				sum += Anrst
						* pow(x_loc, l+n-r)
						* pow(ln_x_loc, s-t);
			}
			Yrt[indexer_Yrt(r, t, ix)]
					= prefactor * sum;
/*if (r==0 and t==0)
	cout << r << "   " << t << "   " << x_loc << "   "
			<< Yrt[indexer_Yrt(r, t, ix)].real() << "   "
			<< Yrt[indexer_Yrt(r, t, ix)].imag() << endl;*/
		}

//if (1) exit (8);

		return;
	}

	void set_ddx_Yrt()
	{
		ddx_Yrt = vector<complex<double> >(rmax*smax*n_x_pts);
		for (int r = 0; r < rmax; ++r)
		for (int t = 0; t < smax; ++t)	//smax is correct
		{
			double h = x_pts[1] - x_pts[0];

			//approximate first derivative at lower endpoint
			ddx_Yrt[indexer_Yrt(r, t, 0)]
				= ( - 3.0 * Yrt[indexer_Yrt(r, t, 0)]
					+ 4.0 * Yrt[indexer_Yrt(r, t, 1)]
					- 1.0 * Yrt[indexer_Yrt(r, t, 2)] )
					/ ( 2.0 * h );

			for (int ix = 1; ix < n_x_pts - 1; ++ix)
			{
				complex<double> fm1
					= Yrt[indexer_Yrt(r, t, ix-1)];
				complex<double> fp1
					= Yrt[indexer_Yrt(r, t, ix+1)];
				ddx_Yrt[indexer_Yrt(r, t, ix)]
					= ( fp1 - fm1 ) / ( 2.0 * h );
			}

			//approximate first derivative at upper endpoint
			ddx_Yrt[indexer_Yrt(r, t, n_x_pts - 1)]
				= ( 3.0 * Yrt[indexer_Yrt(r, t, n_x_pts - 1)]
					- 4.0 * Yrt[indexer_Yrt(r, t, n_x_pts - 2)]
					+ 1.0 * Yrt[indexer_Yrt(r, t, n_x_pts - 3)] )
					/ ( 2.0 * h );
		}

		return;
	}
	//*/


	void get_Yrt(
			vector<complex<double> > * Yrt_in,
			vector<complex<double> > * ddx_Yrt_in,
			vector<double> * x_pts_in,
			int nmax_in, int rmax_in, int l_in)
	{
		nmax = nmax_in;
		rmax = rmax_in;		
		smax = nmax;
		tmax = rmax;
		l = l_in;
		n_x_pts = x_pts_in->size();

		x_pts = *x_pts_in; 

		for (int m = 1; m <= l; ++m)
			sigma_l += 1.0/m;
		//cout << sigma_l << endl;

		B1D = vector<complex<double> >(nmax, 0.0);
		C1D = vector<complex<double> >(nmax, 0.0);

		A3D = vector<complex<double> >(rmax*nmax*smax, 0.0);
		B3D = vector<complex<double> >(rmax*nmax*smax, 0.0);
		C3D = vector<complex<double> >(rmax*nmax*smax, 0.0);

		//cout << "Start." << endl;
		initialize_A3D_starting_values();

		//cout << "Next." << endl;
		for (int n = 2; n < nmax; ++n)
		{
			cerr << "n = " << n << endl;
			set_B1D(n);
			set_C1D(n);

			for (int t = 0; t < rmax; ++t)
			for (int s = 0; s < smax; ++s)
			{
				set_B3D(n, t, s);
				set_C3D(n, t, s);
			}

			for (int r = 0; r < rmax; ++r)
			for (int s = 0; s < smax; ++s)
				set_A3D(n, r, s);
		}
		//cout << "Here." << endl;

/*
		//checks the X_n(y) functions
		for (int n = 0; n < nmax; ++n)
		for (int y = 0; y < 101; ++y)
		{
			complex<double> sum = 0.0;
			double y0 = 1.0 + 0.01*y;
			for (int r = 0; r < rmax; ++r)
			for (int s = 0; s < smax; ++s)
				sum += A3D.at(indexer_A3D(n, r, s))
						* pow(y0, l+n-r) * pow(log(y0), s);
			if ( n > 0 and abs(y0-1.0) < 1.e-30 ) 
				sum = 0.0;
			
			cout << n << "   " << y0 << "   "
				<< setprecision(15)
				<< sum.real() << "   " << sum.imag()
				<< endl;
		}
		cout << endl;

if (1) exit (8);
*/
///*
		//checks the X_n(y) functions
		for (int y = 0; y < 91; ++y)
		{
			complex<double> sum = 0.0;
			double y0 = 1.0 + 0.1*y;
			double xs = 1.0;		//for simplicity, for checking purposes
			double x = y0 * xs;
			vector<double> inner_sums_re, inner_sums_im;
			for (int n = 0; n < nmax; ++n)
			{
				complex<double> inner_sum = 0.0;
				for (int r = 0; r < rmax; ++r)
				for (int s = 0; s < smax; ++s)
					sum += double(r==0)*A3D.at(indexer_A3D(n, r, s))
							* pow(y0, l+n-r)
							* pow(log(y0), s)
							* pow(x, n);
				sum += inner_sum;
				if (n%2==0)
					inner_sums_re.push_back(inner_sum.real());
				else
					inner_sums_im.push_back(inner_sum.imag());
			}

			//this accelerate real and imaginary parts of sum separately
			double sum_accel_re = 0.0, sum_accel_im = 0.0, err = 0.0;
			if (USE_SERIES_ACCELERATION and false)
			{
				compute_levin_sum(inner_sums_re, &sum_accel_re, &err);
				compute_levin_sum(inner_sums_im, &sum_accel_im, &err);

				sum = complex<double> (sum_accel_re, sum_accel_im);
			}
			
			cout << y0 << "   "
				<< setprecision(15)
				<< sum.real() << "   " << sum.imag()
				<< endl;
		}
		cout << endl;

if (1) exit (8);
//*/

		/*for (int n = 0; n < nmax; ++n)
		for (int r = 0; r < rmax; ++r)
		for (int s = 0; s < smax; ++s)
		{
			complex<double> A_nrs = A3D.at(indexer_A3D(n, r, s));
			for (int t = 0; t <= s; ++t)	//smax is correct
				A4D.at( indexer_A4D(n, r, s, t) )
					= A_nrs * pow(-1.0,t) * gsl_sf_choose(s, t);
		}*/

		set_Yrt();
		
		set_ddx_Yrt();

		//finally, be sure to copy over results!
		for (int iYrt = 0; iYrt < Yrt.size(); ++iYrt)
			Yrt_in->at(iYrt) = Yrt.at(iYrt);

		for (int iYrt = 0; iYrt < ddx_Yrt.size(); ++iYrt)
			ddx_Yrt_in->at(iYrt) = ddx_Yrt.at(iYrt);

		return;
	}
}

