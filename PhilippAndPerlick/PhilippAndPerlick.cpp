#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

using namespace std;

#include "RecurrenceRelation.h"

void linspace(vector<double> & x, double a, double b);

const complex<double> i(0.0, 1.0);
const double epsilon = 1.e-10;
bool check_for_convergence = true;

struct rec_rel_params
{
	complex<double> Omega1;
	complex<double> Omega2;
	complex<double> Omega3;
	complex<double> Delta1;
	complex<double> Delta2;
};

//define the coefficient functions
inline complex<double> alpha(int n, void * p);
inline complex<double> beta(int n, void * p);
inline complex<double> gamma(int n, void * p);
inline complex<double> a(int n, void * p);
inline complex<double> b(int n, void * p);

int main (int argc, char* argv[])
{
	const long int N = 5000;
	const long int kmax = 5000;

	//set initial parameters
	const complex<double> xi0 = 1.0;			//arbitrary, apparently
	complex<double> l = 0.0;						//ang. mom. of emitted wave
	complex<double> s = 0.0;								//spin of scalar field s==0
	//the frequency of the emitted wave
	//complex<double> omega = 10.0;
	complex<double> omega = complex<double>(atof(argv[1]));

	//set indices to specify given solution
	//complex<double> nu = 2.0*i*omega;
	//complex<double> mu1 = 2.0*i*omega;
	//complex<double> mu2 = -(s+1.0);
	complex<double> nu = 0.0;
	complex<double> mu1 = 0.0;
	complex<double> mu2 = -(s+1.0);

	//set r grid at which to compute function
	const int n_r_pts = 2;
	vector<double> r_pts(n_r_pts);
	linspace(r_pts, 200.0, 5000.0);

	//set coefficients in recurrence relation
	complex<double> Omega1 = 1.0 - 2.0*i*omega + 2.0*mu1;
	complex<double> Omega2 = 2.0 * (s + nu + mu2 - mu1 - 1.0);
	complex<double> Omega3 = 1.0 - 2.0*s - 2.0*mu2;
	complex<double> Delta1 = s * (s + 1.0) - l * (l + 1.0) + 2.0 * (s + mu2)*(mu1 - nu)
								- nu + mu1 + mu2;
	complex<double> Delta2 = mu2 * (2.0 * s + mu2);

	struct rec_rel_params params = { Omega1, Omega2, Omega3, Delta1, Delta2 };

//cout << "CHECK: " << omega << "   " << nu << "   " << mu1 << "   " << mu2 << "   " << Omega1 << "   " << Omega2 << "   " << Omega3 << "   " << Delta1 << "   " << Delta2 << endl;

	//initialize parameters for recurrence relation solution
	RecurrenceRelation::threshold = epsilon;
	RecurrenceRelation::N = N;
	RecurrenceRelation::kmax = kmax;

	//loop over r grid, get minimal solution to recurrence relation at each point
	for (int ir = 0; ir < n_r_pts; ++ir)
	{
		double r = r_pts[ir];
		complex<double> prefactor = exp(nu*(r-1.0))*pow(r-1.0, mu1)*pow(r,mu2+s+1.0);
		
		//vectors to hold ratios of consecutive terms in minimal
		//solution and the minimal solution itself, respectively
		//minimal solution computed in range f_0..f_N
		vector<complex<double> > rn(N+1), results(N+1);

		//set the last ratio first, along with initial term in series
		results[0] = xi0;
		complex<double> xi1 = -beta(0, &params)*xi0/alpha(0, &params);
		results[1] = xi1;
		results[2] = -(beta(1, &params)*xi1+gamma(1, &params)*xi0)/alpha(1, &params);
		rn[N] = RecurrenceRelation::get_rN(&a, &b, &params);

		//now get the actual solution (this algorithm described by Gautschi (1967))
		RecurrenceRelation::solve_recurrence_relation(&rn, &results, &a, &b, &params);

		bool converged = false;
		complex<double> sum = 0.0;
		complex<double> xik = 0.0;
		double r_coefficient = 1.0;
		for (int k = 0; k < results.size(); ++k)
		{
			xik = results[k];
			sum += xik * r_coefficient;
			if ( check_for_convergence and abs(xik * r_coefficient) / abs(sum) < epsilon )
			{
				converged = true;
				break;	//converged; terminate loop
			}
			r_coefficient *= (r-1.0)/r;
			if (k<=5 and ir == n_r_pts - 1)
				cout << "TEST: " << k << "   " << setprecision(8)
						<< "alpha_" << k << " xi_" << k+1 << " + beta_"
						<< k << " xi_" << k << " + gamma_"
						<< k << " xi_" << k-1 << " = "
						<< results[k+1]*alpha(k, &params)
							+ results[k]*beta(k, &params)
							+ results[k-1]*gamma(k, &params)
						<< endl << "\t -->rn[" << k << "] = " << rn[k] << "   " << results[k+1]/results[k]
						 << "   " << -beta(0, &params)/alpha(0, &params) << endl;
		}

		sum *= prefactor;
		if (check_for_convergence and not converged)
			cerr << "Did not converge at r = " << r
					<< "!!  effective epsilon = " << 2.0*abs(xik * r_coefficient) / abs(sum)
					<< "; " << abs(xik * r_coefficient) << "   " << abs(sum) << endl;
		if (ir == n_r_pts - 1)
			cout << setw(16) << setprecision(12)
					<< omega.real() << "   " << r << "   "
					<< sum.real() << "   " << sum.imag()
					<< "   " << 1.0/( abs(sum) * abs(sum) ) << endl;
	}
	
	return 0;
}

void linspace(vector<double> & x, double a, double b)
{
//returns vector x of n linearly spaced points, with a and b as endpoints
	int n = x.size();
	double Del_x = (b-a)/(double)(n-1);

	//assume x has length n already
	for (int i = 0; i < n; i++)
		x[i] = a + Del_x * (double)i;

	return;
}

//define the coefficient functions
inline complex<double> alpha(int k, void * p)
{
	struct rec_rel_params * params = (struct rec_rel_params *)p;
	complex<double> Omega1 = (params->Omega1);
	double dk = k;
	//complex<double> k_0_result = Omega1 + 1.0;
	complex<double> k_0_result = Omega1;
	complex<double> k_pos_result = 1.0 + (Omega1 + 1.0)/dk + Omega1 / (dk*dk);

	if (k==0)
		return k_0_result;
	else
		return k_pos_result;

	return ( double(k==0)*k_0_result + double(k>0)*k_pos_result );
}

inline complex<double> beta(int k, void * p)
{
	struct rec_rel_params * params = (struct rec_rel_params *)p;
	complex<double> Omega2 = (params->Omega2);
	complex<double> Delta1 = (params->Delta1);
	double dk = k;
	complex<double> k_0_result = Delta1;
	complex<double> k_pos_result = -2.0 + (Omega2 + 2.0)/dk + Delta1 / (dk*dk);

	if (k==0)
		return k_0_result;
	else
		return k_pos_result;

	return ( double(k==0)*k_0_result + double(k>0)*k_pos_result );
}

inline complex<double> gamma(int k, void * p)
{
	struct rec_rel_params * params = (struct rec_rel_params *)p;
	complex<double> Omega3 = (params->Omega3);
	complex<double> Delta2 = (params->Delta2);
	double dk = k;
	complex<double> k_0_result = 0.0;
	//complex<double> k_pos_result = 1.0 - (Omega3 - 3.0)/dk + (Delta2 - Omega3 + 2.0) / (dk*dk);
	complex<double> k_pos_result = 1.0 + (Omega3 - 3.0)/dk + (Delta2 - Omega3 + 2.0) / (dk*dk);

	if (k==0)
		return k_0_result;
	else
		return k_pos_result;

	return ( double(k==0)*k_0_result + double(k>0)*k_pos_result );
}

inline complex<double> a(int n, void * p)
{
	return( beta(n, p)/alpha(n, p) );
}

inline complex<double> b(int n, void * p)
{
	return( gamma(n, p)/alpha(n, p) );
}

