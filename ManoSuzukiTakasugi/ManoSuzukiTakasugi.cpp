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
const int l = 2, s = 2;

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
complex<double> discriminant(void * p);
complex<double> evaluate_alpha0R1_continued_fraction(
					complex<double>(*alpha_in)(int,void*),
					complex<double>(*beta_in)(int,void*),
					complex<double>(*gamma_in)(int,void*),
					void * p );
complex<double> evaluate_gamma0Lm1_continued_fraction(
					complex<double>(*alpha_in)(int,void*),
					complex<double>(*beta_in)(int,void*),
					complex<double>(*gamma_in)(int,void*),
					void * p );


int main (int argc, char* argv[])
{
	for (int inu = 0; inu <= 1000; inu++)
	{
		double nu = 1.9 + 0.0001*inu;
		struct function_params params = { nu, 1.0/sqrt(10.0) };

		complex<double> result = discriminant(&params);

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

	return ( numerator / (denominator+0.e-100) );
}

inline complex<double> beta(int n, void * p)
{
	struct function_params * params = (struct function_params *)p;
	double nu = (params->nu);
	double e = (params->epsilon);
	double e2 = e*e;

	return ( (n+nu)*(n+nu+1.0+0.e-100) - l*(l+1.0) + 2.0*e2
				+ e2*(s*s+e2)/((n+nu)*(n+nu+1.0)+0.e-100) );
}

inline complex<double> gamma(int n, void * p)
{
	struct function_params * params = (struct function_params *)p;
	double nu = (params->nu);
	double e = (params->epsilon);
	complex<double> f3 = (n+nu-i*e)*(n+nu-i*e);

	complex<double> numerator = i*e*(n+nu+2.0+i*e)*(n+nu+i*e)*(f3-double(s*s));
	complex<double> denominator = (n+nu-2.0-i*e)*(n+nu)*(2.0*n+2.0*nu-1.0);

	return ( numerator / (denominator+0.e-100) );
}

complex<double> discriminant(void * p)
{
	return (
			beta(0, p)
			+ evaluate_alpha0R1_continued_fraction( &alpha, &beta, &gamma, p )
			+ evaluate_gamma0Lm1_continued_fraction( &alpha, &beta, &gamma, p )
			);
}

//uses modified Lentz's method (cf. Numerical Recipes in C++)
complex<double> evaluate_alpha0R1_continued_fraction(
					complex<double>(*alpha_in)(int,void*),
					complex<double>(*beta_in)(int,void*),
					complex<double>(*gamma_in)(int,void*),
					void * p )
{
	const int jmax = 100;
	const double tiny = tolerance*1.e-20;
	complex<double> result = tiny;
	complex<double> C = result, D = 0.0;
	const int start_index = 1;

	for (int j = 1; j <= jmax; j++)
	{
		//cerr << "BEFORE: " << j << "   " << C << "   " << D << "   " << result << endl;
		//be careful with indexing
		complex<double> a = -alpha(start_index + j - 2, p)
							* gamma(start_index + j - 1, p);
		complex<double> b = beta(start_index + j - 1, p);

		D = b + a * D;
		if (abs(D) < tiny)
			D = tiny;

		C = b + a / C;
		if (abs(C) < tiny)
			C = tiny;

		D = 1.0 / D;

		complex<double> Delta = C * D;
		result *= Delta;
		if ( abs(Delta - 1.0) < tolerance )
			break;
		else if (j == jmax)
			cerr << "Warning: modified Lentz's method did not converge!" << endl
					<< "alpha0R1(): abs(Delta - 1.0) = " << abs(Delta - 1.0) 
					<< " >= tolerance = " << tolerance << endl;
		//cerr << "AFTER: " << j << "   " << C << "   " << D << "   "
		//		<< Delta << "   " << Delta - 1.0 << "   " << abs(Delta - 1.0) << "   "
		//		<< result << endl;
	}

	return (result);
}


//uses modified Lentz's method (cf. Numerical Recipes in C++)
complex<double> evaluate_gamma0Lm1_continued_fraction(
					complex<double>(*alpha_in)(int,void*),
					complex<double>(*beta_in)(int,void*),
					complex<double>(*gamma_in)(int,void*),
					void * p )
{
	const int jmax = 100;
	const double tiny = tolerance*1.e-20;
	complex<double> result = tiny;
	complex<double> C = result, D = 0.0;
	const int start_index = -1;

	for (int j = 1; j <= jmax; j++)
	{
		//cerr << "BEFORE: " << j << "   " << C << "   " << D << "   " << result << endl;
		//be careful with indexing
		complex<double> a = -alpha(start_index - j + 1, p)
							* gamma(start_index - j + 2, p);
		complex<double> b = beta(start_index - j + 1, p);

		D = b + a * D;
		//cerr << "check: " << j << "   " << C << "   " << D << "   "
		//		<< alpha(start_index - j + 1, p) << "   " << gamma(start_index - j + 2, p) << "   "
		//		<< result << endl;
		if (abs(D) < tiny)
			D = tiny;
		//cerr << "check2: " << j << "   " << C << "   " << D << "   " << result << endl;

		C = b + a / C;
		//cerr << "check3: " << j << "   " << C << "   " << D << "   " << result << endl;
		if (abs(C) < tiny)
			C = tiny;
		//cerr << "check4: " << j << "   " << C << "   " << D << "   " << result << endl;

		D = 1.0 / D;
		//cerr << "check5: " << j << "   " << C << "   " << D << "   " << result << endl;

		complex<double> Delta = C * D;
		result *= Delta;
		if ( abs(Delta - 1.0) < tolerance )
			break;
		else if (j == jmax)
			cerr << "Warning: modified Lentz's method did not converge!" << endl
					<< "alpha0Lm1(): abs(Delta - 1.0) = " << abs(Delta - 1.0) 
					<< " >= tolerance = " << tolerance << endl;
		//cerr << "AFTER: " << j << "   " << C << "   " << D << "   "
		//		<< Delta << "   " << Delta - 1.0 << "   " << abs(Delta - 1.0) << "   "
		//		<< result << endl;
	}

	return (result);
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

