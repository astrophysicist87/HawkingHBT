#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

using namespace std;

#include "RecurrenceRelation.h"

//some parameters
const complex<double> i(0.0,1.0);
const double threshold = 1.e-10;
const int N = 10, kmax = 10;
//const double BesselJ_0_at_z = 0.76519768655796655144971752610266322090927428975533;
const complex<double> BesselJ_0_at_z = -0.5936889765458604437 - 4.8838405750713968930*i;

struct my_f_params { complex<double> z; };

//define the coefficient functions
inline complex<double> alpha(int n, void * p)
{
	struct my_f_params * params = (struct my_f_params *)p;
	complex<double> z = (params->z);
	return (z);
}

inline complex<double> beta(int n, void * p)
{
	struct my_f_params * params = (struct my_f_params *)p;
	complex<double> z = (params->z);
	return (-2.0*n);
}

inline complex<double> gamma(int n, void * p)
{
	struct my_f_params * params = (struct my_f_params *)p;
	complex<double> z = (params->z);
	return (z);
}

inline complex<double> a(int n, void * p)
{
	return( beta(n, p)/alpha(n, p) );
}

inline complex<double> b(int n, void * p)
{
	return( gamma(n, p)/alpha(n, p) );
}

//the main driver function
int main(void)
{
	//N+1 is the number of terms of minimal solution
	//to be computed, including 0th term
	vector<complex<double> > rn(N+1), results(N+1);
	results[0] = BesselJ_0_at_z;
	complex<double> chosen_z = 2.0 + M_PI*i;

	struct my_f_params params = { chosen_z };

	RecurrenceRelation::threshold = threshold;
	RecurrenceRelation::N = N;
	RecurrenceRelation::kmax = kmax;

	rn[N] = RecurrenceRelation::get_rN(&a, &b, &params);

	RecurrenceRelation::solve_recurrence_relation(&rn, &results, &a, &b, &params);

	for (int iN = 0; iN < results.size(); ++iN)
		cout << setw(14) << setprecision(12)
				<< "J_" << iN << "(" << chosen_z << ") = " << results[iN] << endl;

	return 0;
}
