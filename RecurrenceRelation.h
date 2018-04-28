#ifndef RECURRENCE_RELATION_H
#define RECURRENCE_RELATION_H

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

using namespace std;

namespace RecurrenceRelation
{
	//some parameters
	const complex<double> i(0.0,1.0);
	double threshold = 1.e-1;
	int N = -10, kmax = -10;

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
	}}

#endif
