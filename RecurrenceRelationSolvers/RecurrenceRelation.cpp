#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

using namespace std;

//some parameters
const complex<double> i(0.0,1.0);
const double z = 1.0, threshold = 1.e-10;
const int N = 10, kmax = 5;	//probably good enough
const double BesselJ_0_at_1 = 0.76519768655796655144971752610266322090927428975533;

//define the coefficient functions
inline double a(int n)
{
	return (-2.0 * n / z);
}

inline double b(int n)
{
	return (1.0);
}

//the main driver function
int main(void)
{
	double u = 1.0, v = -b(N+1)/a(N+1);
	double w = v;
	double relative_increment = 1.0;

	int k = 1;
	do
	{
		u = 1.0 / (1.0 - u*b(N+k+1)/(a(N+k)*a(N+k+1)));
		v *= u - 1.0;
		relative_increment = v/w;
		w += v;
		//cout << "r_" << k << " = " << w << endl;
		k++;
	} while (abs(relative_increment) >= threshold and k <= kmax);

	vector<double> rn(N+1), results(N+1);
	results[0] = BesselJ_0_at_1;

	rn[N] = w;
	for (int iN = N; iN >= 1; --iN)
		rn[iN-1] = -b(iN)/(a(iN)+rn[iN]);
	
	for (int iN = 1; iN < results.size(); ++iN)
		results[iN] = rn[iN-1]*results[iN-1];

	for (int iN = 0; iN < results.size(); ++iN)
		cout << setw(14) << setprecision(12) << "J_" << iN << "(" << z << ") = " << results[iN] << "; r[" << iN << "] = " << rn[iN] << endl;

	return 0 ;
}
