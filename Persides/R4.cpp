#include <stdio.h>
#include <iostream>
#include <complex>
#include <vector>

#include <gsl/gsl_sf_gamma.h>

using namespace std;

vector<complex<double> > B1D, C1D;
vector<complex<double> > A3D, B3D, C3D, D3D;
vector<complex<double> > A4D;

void initialize_A3D_starting_values();

int main(void)
{
	const int rmax = 50, nmax = 50, smax = 50, tmax = 50;

	double sigma_l = 1.0;
	if (l >= 1)
	{
		sigma_l = 0.0;
		for (int m = 1; m <= l; ++m)
			sigma_l += 1.0/m;
	}

	A3D = vector<complex<double> >(rmax*nmax*smax, 0.0);

	initialize_A3D_starting_values();
	
	for (int n = 2; n <= nmax; ++n)
	{
		set_B1D(n);
		set_C1D(n);

		set_B3D(n);
		set_C3D(n);

		set_A3D(n);
	}

	return 0 ;
}

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

void initialize_A3D_starting_values()
{
	//set n==0,s==0
	for (int r = 0; r < rmax; ++r)
		A3D[indexer_A3D(0, r, 0)]
			= pow(-1.0, l)*pp(r);

	//set n==1,s==0
	for (int r = 0; r < rmax; ++r)
	{
		double small_sum = 0.0;
		for (int m = 0; m <= r-2; ++m)
			small_sum += pp(r-m-2)/(m+1.0);
		
		A3D[indexer_A3D(1, r, 0)]
			= (1.0-2.0*sigma_l)*i*pow(-1.0,l)*pp(r-1)
				+ 2.0*i*pow(-1.0,l)*qp(r-2*l-2)
				- i*pow(-1.0,l)*small_sum;
	}

	//set n==1,s==1
	for (int r = 0; r < rmax; ++r)
		A3D[indexer_A3D(1, r, 1)]
			= i*pow(-1.0, l+1)*pp(r-1);

	return;
}

void set_B1D(int n)
{
	double sum = 0.0;
	for (int u = 0; u < rmax; ++u)
	for (int t = 0; t <= u+l; ++t)	//double check
	for (int v = 0; v < smax; ++v)
		sum += (l-t+u+1.0)*pp(t-u)*A3D[indexer_A3D(n-1,u,v)]*D3D(2*l+n-1,v,0);

	B1D[n] = sum;

	return;
}

void set_C1D(int n)
{
	double sum = 0.0;
	for (int t = 0; t < tmax; ++t) //double check
	for (int u = 0; u < rmax; ++u)
	for (int v = 0; v < smax; ++v)
		sum += (-l-t+u)*qp(t-u)*A3D[indexer_A3D(n-1,u,v)]*D3D(n-1-t,v,0);

	C1D[n] = sum;

	return;
}

void set_B3D(int n, int t, int s)
{
	double sum = 0.0;
	for (int u = 0; u < rmax; ++u)
	for (int v = 0; v < smax; ++v)
		sum += (l-t+u+1.0)*pp(t-u)*A3D[indexer_A3D(n-1,u,v)]*D3D(2*l+n-t,v,s);

	B3D[indexer_B3D(n, t, s)] = sum - B1D[n]*double(s==0)*double(2*l+n+1==t);

	return;
}

void set_C3D(int n, int t, int s)
{
	double sum = 0.0;
	for (int u = 0; u < rmax; ++u)
	for (int v = 0; v < smax; ++v)
		sum += (-l-t+u)*qp(t-u)*A3D[indexer_A3D(n-1,u,v)]*D3D(n-1-t,v,s);

	C3D[indexer_C3D(n, t, s)] = sum - C1D[n]*double(s==0)*double(n==t);

	return;
}

void set_A3D(int n, int r, int s)
{
	double sum = 0.0;
	for (int t = 0; t <= r; ++t) //double check
		sum += qp(r-t)*B3D[indexer_B3D(n, t, s)]
				- pp(r-t)*C3D[indexer_C3D(n, t, s)];

	A3D[indexer_A3D(n,r,s)] = 4.0*i*sum;
}

void set_Yrt()
{
	double prefactor = pow(-1.0,l)*Al(l);
	for (int ix = 0; ix < n_x_pts; ++ix)
	for (int r = 0; r < rmax; ++r)
	for (int t = 0; t < tmax; ++t)
	{
		if (r<t)
			continue;
		double x_loc = x_pts[ix];
		double ln_x_loc = log(x_loc);
		complex<double> sum = 0.0;
		for (int n = 0; n < nmax; ++n)
		for (int s = 0; s < smax; ++s)
		{
			if (t>s)
				continue;
			double Anrst = pow(-1.0,t)*gsl_sf_choose(s, t)*A3D[indexer_A3D(n,r,s)];
			sum += Anrst*pow(x_loc, l+n-r)*pow(ln_x_loc, s-t);
		}
		Yrt[indexer_Yrt(r, t, ix)] = prefactor * sum;
	}

	return;
}

void set_xinm(int n, int m)
{
	if (n<m)
		return;
	else
	{
		for (int ix = 0; ix < n_x_pts; ++ix)
		{
			complex<double> sum = 0.0;
			for (int s = 0; s < smax; ++s)
			{
			
			}
			xinm[indexer_xinm(t, m, ix)] = conj(sum);
		}
	}
	return;
}
