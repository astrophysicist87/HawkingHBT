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
	const int rmax = 50, nmax = 50, smax = 50;

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
		set_D3D(n);

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
