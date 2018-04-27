#include <stdio.h>
#include <iostream>
#include <complex>

using namespace std;

const complex<double> i(0.0,1.0);

void run_asymptotics_v1();

int main(void)
{
	run_asymptotics_v1();

	return 0 ;
}

void run_asymptotics_v1()
{
	int l = 0;
	int nmax = 100;
	double x = 5.0, xs = 1.0;
	complex<double> tau0 = 1.0;
	complex<double> tau1 = 0.5*i*double(l)*(l+1.0)*tau0;

	complex<double> tau_nm2 = tau0;
	complex<double> tau_nm1 = tau1;

	complex<double> sum = tau0/x;
	cout << 0 << "   " << sum.real() << "   " << sum.imag() << endl;

	sum += tau1/(x*x);
	cout << 1 << "   " << sum.real() << "   " << sum.imag() << endl;

	double x_power = x*x;
	
	for (int n = 2; n <= nmax; ++n)
	{
		complex<double> tau_n = i*( (l+n)*(l-n+1.0)*tau_nm1 + (n-1.0)*(n-1.0)*xs*tau_nm2 )/(2.0*n);
		x_power *= x;

		sum += tau_n / x_power;

		tau_nm2 = tau_nm1;
		tau_nm1 = tau_n;

		cout << n << "   " << sum.real() << "   " << sum.imag() << endl;
	}

}
