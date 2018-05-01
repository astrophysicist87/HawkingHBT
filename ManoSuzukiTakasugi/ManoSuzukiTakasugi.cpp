#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

#include <gsl/gsl_sf_gamma.h>

using namespace std;

#include "RecurrenceRelation.h"

void linspace(vector<double> & x, double a, double b);

const complex<double> i(0.0, 1.0);
const double tolerance = 1.e-10;
bool check_for_convergence = true;
const int l = 0, s = 0;

struct function_params
{
	double nu;
	double epsilon;
};

//define the coefficient functions
inline complex<double> alpha(int n, void * p);
inline complex<double> beta(int n, void * p);
inline complex<double> gamma(int n, void * p);
inline complex<double> a(int n, void * p);
inline complex<double> b(int n, void * p);
complex<double> get_discriminant(void * p);
complex<double> get_R1(complex<double>(*a_in)(int,void*),
					complex<double>(*b_in)(int,void*),
					void * p );
complex<double> get_L0(complex<double>(*a_in)(int,void*),
					complex<double>(*b_in)(int,void*),
					void * p );


int main (int argc, char* argv[])
{
	for (int inu = 0; inu <= 1000; inu++)
	{
		double nu = -0.5 + 0.001*inu;
		struct function_params params = { nu, 0.1 };

		complex<double> result = get_discriminant(&params);

		cout << nu << "   " << result.real() << "   " << result.imag() << endl;
	}
	
	return 0;
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

/*
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
*/

//define the coefficient functions
inline complex<double> alpha(int n, void * p)
{
	struct function_params * params = (struct function_params *)p;
	double nu = (params->nu);
	double e = (params->epsilon);
	complex<double> f3 = (n+nu+1.0+i*e)*(n+nu+1.0+i*e);

	complex<double> numerator = -i*e*(n+nu+1.0-i*e)*(n+nu-1.0-i*e)*(f3-double(s*s));
	complex<double> denominator = (n+nu+3.0+i*e)*(n+nu+1.0)*(2.0*n+2.0*nu+3.0);

	return ( numerator / denominator );
}

inline complex<double> beta(int n, void * p)
{
	struct function_params * params = (struct function_params *)p;
	double nu = (params->nu);
	double e = (params->epsilon);
	double e2 = e*e;

	return ( (n+nu)*(n+nu+1.0) - l*(l+1.0) + 2.0*e2
				+ e2*(s*s+e2)/((n+nu)*(n+nu+1.0)) );
}

inline complex<double> gamma(int n, void * p)
{
	struct function_params * params = (struct function_params *)p;
	double nu = (params->nu);
	double e = (params->epsilon);
	complex<double> f3 = (n+nu-i*e)*(n+nu-i*e);

	complex<double> numerator = i*e*(n+nu+2.0+i*e)*(n+nu+i*e)*(f3-double(s*s));
	complex<double> denominator = (n+nu-2.0-i*e)*(n+nu)*(2.0*n+2.0*nu-1.0);

	return ( numerator / denominator );
}

inline complex<double> a(int n, void * p)
{
	return( beta(n, p)/alpha(n, p) );
}

inline complex<double> b(int n, void * p)
{
	return( gamma(n, p)/alpha(n, p) );
}



complex<double> get_discriminant(void * p)
{
	complex<double> R1 = get_R1(&a, &b, p);
	complex<double> L0 = get_L0(&a, &b, p);

	return (R1*L0 - 1.0);
}

complex<double> get_R1(complex<double>(*a_in)(int,void*),
					complex<double>(*b_in)(int,void*),
					void * p )
{
	double threshold = 1.e-10;
	const int N = 0, kmax = 10;
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

complex<double> get_L0(complex<double>(*a_in)(int,void*),
					complex<double>(*b_in)(int,void*),
					void * p )
{
	double threshold = 1.e-10;
	const int N = -1, kmax = 10;
	complex<double> u = 1.0, v = -b_in(N-1,p)/a_in(N-1,p);
	complex<double> w = v;
	complex<double> relative_increment = 1.0;

	int k = 1;
	do
	{
		u = 1.0 / (1.0 - u*b_in(N-k-1,p)/(a_in(N-k,p)*a_in(N-k-1,p)));
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



/*
complex<double> pHypGeo(double x, void * p)
{
	struct function_params * params = (struct function_params *)p;
	double n = (params->n);
	double nu = (params->nu);
	double epsilon = (params->epsilon);

	gsl_sf_result G1_lnr, G1_arg, G2_lnr, G2_arg, G3_lnr, G3_arg;

	gsl_sf_lngamma_complex_e (n+nu-1.0, -epsilon, &G1_lnr, &G1_arg);
	gsl_sf_lngamma_complex_e (-n-nu-2.0, -epsilon, &G2_lnr, &G2_arg);
	gsl_sf_lngamma_complex_e (1.0, -2.0*epsilon, &G3_lnr, &G3_arg);
	complex<double> prefactor = exp(G1_lnr + G2_lnr - G3_lnr + i*(G1_arg + G2_arg - G3_arg));

	return ( prefactor * Hypergeometric2F1(n+nu-1.0-i*epsilon, -n-nu-2.0-i*epsilon, 1.0-i*2.0*epsilon, x) );
}
*/

