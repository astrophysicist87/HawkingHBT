#include <stdio.h>
#include <iostream>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sum.h>

using namespace std;

void compute_levin_sum(vector<double> computed_terms, double * sum_accel, double * err);

int main (void)
{
	const int nterms = 20;
	const int starting_term = 5;

	vector<double> computed_terms(nterms, 0.0);
	for (int n = 0; n < nterms; n++)
	{
		double np1 = starting_term + n + 1.0;
		computed_terms[n] = 1.0 / (np1*np1);
	}

	double sum_accel, err;
	compute_levin_sum(computed_terms, &sum_accel, &err);
	cout << sum_accel << "   " << err << endl;

	return (0);
}

void compute_levin_sum(vector<double> computed_terms, double * sum_accel, double * err)
{
	int N = computed_terms.size();
	double * t = new double[N];
	double sum = 0;

	gsl_sum_levin_u_workspace * w = gsl_sum_levin_u_alloc (N);

	//const double zeta_2 = M_PI * M_PI / 6.0;

	for (int n = 0; n < N; n++)
	{
		t[n] = computed_terms[n];
		sum += t[n];
	}

	gsl_sum_levin_u_accel (t, N, w, sum_accel, err);

	//printf ("term-by-term sum = % .16f using %d terms\n", sum, N);

	//printf ("term-by-term sum = % .16f using %zu terms\n", w->sum_plain, w->terms_used);

	//printf ("exact value      = % .16f\n", zeta_2);
	//printf ("accelerated sum  = % .16f using %zu terms\n", sum_accel, w->terms_used);

	//printf ("estimated error  = % .16f\n", err);
	//printf ("actual error     = % .16f\n", sum_accel - zeta_2);

	delete [] t;

	gsl_sum_levin_u_free (w);

	return;
}
