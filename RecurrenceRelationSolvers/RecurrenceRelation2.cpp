#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

using namespace std;

//some parameters
const complex<double> i(0.0,1.0);
const double threshold = 1.e-10;
const int N = 10, kmax = 10;	//probably good enough
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

//function prototypes
complex<double> get_rN(
			complex<double>(*a_in)(int,void*),
			complex<double>(*b_in)(int,void*),
			void * p );
void solve_recurrence_relation(
			vector<complex<double> > * rn_ptr,
			vector<complex<double> > * results_ptr,
			complex<double>(*a_in)(int,void*),
			complex<double>(*b_in)(int,void*),
			void * p );

//the main driver function
int main(void)
{
	//N+1 is the number of terms of minimal solution
	//to be computed, including 0th term
	vector<complex<double> > rn(N+1), results(N+1);
	results[0] = BesselJ_0_at_z;
	complex<double> chosen_z = 2.0 + M_PI*i;

	struct my_f_params params = { chosen_z };

	rn[N] = get_rN(&a, &b, &params);

	solve_recurrence_relation(&rn, &results, &a, &b, &params);

	for (int iN = 0; iN < results.size(); ++iN)
		cout << setw(14) << setprecision(12)
				<< "J_" << iN << "(" << chosen_z << ") = " << results[iN] << endl;

	return 0 ;
}

complex<double> get_rN(complex<double>(*a_in)(int,void*),
				complex<double>(*b_in)(int,void*),
				void * p )
{
	complex<double> u = 1.0, v = -b_in(N+1,p)/a_in(N+1,p);
	complex<double> w = v;
	complex<double> relative_increment = 1.0;

	int k = 1;
	do
	{
		u = 1.0 / (1.0 - u*b_in(N+k+1,p)/(a_in(N+k,p)*a_in(N+k+1,p)));
		v *= u - 1.0;
		relative_increment = v/w;
		w += v;
		//cout << "r_" << k << " = " << w << endl;
		k++;
	} while (abs(relative_increment) >= threshold and k <= kmax);
	if (abs(relative_increment) >= threshold)
		cerr << "Warning in get_rN(): did not converge to desired threshold!" << endl
				<< "\t abs(relative_increment) == " << abs(relative_increment) << endl
				<< "\t threshold == " << threshold
				<< endl;

	return (w);
}

void solve_recurrence_relation(
			vector<complex<double> > * rn_ptr,
			vector<complex<double> > * results_ptr,
			complex<double>(*a_in)(int,void*),
			complex<double>(*b_in)(int,void*),
			void * p )
{
	for (int iN = N; iN >= 1; --iN)
		rn_ptr->at(iN-1) = -b_in(iN,p)/(a_in(iN,p)+rn_ptr->at(iN));
	
	for (int iN = 1; iN < results_ptr->size(); ++iN)
		results_ptr->at(iN) = rn_ptr->at(iN-1)*results_ptr->at(iN-1);

	return;
}
