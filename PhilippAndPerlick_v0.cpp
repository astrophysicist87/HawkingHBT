#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

using namespace std;

void linspace(vector<double> & x, double a, double b);

const complex<double> i(0.0, 1.0);
const double epsilon = 1.e-12;
bool check_for_convergence = true;

int main (void)
{
	const long int kmax = 10000;
	const complex<double> xi0 = 1.0;			//arbitrary, apparently

	complex<double> l = 30.0;						//ang. mom. of emitted wave
	complex<double> s = 0.0;								//spin of scalar field s==0
	complex<double> omega = 30.0;						//the frequency of the emitted wave
												//- not sure what a reasonable value is

	complex<double> nu = 2.0*i*omega;
	complex<double> mu1 = 2.0*i*omega;
	complex<double> mu2 = -(s+1.0);

	const int n_r_pts = 1000;
	vector<double> r_pts(n_r_pts);
	linspace(r_pts, 1.5, 2.0);

	for (int ir = 0; ir < n_r_pts; ++ir)
	{
		double r = r_pts[ir];
		complex<double> prefactor = exp(nu*(r-1.0))*pow(r-1.0, mu1)*pow(r,mu2+s+1.0);
		
		complex<double> Omega1 = 1.0 - 2.0*i*omega + 2.0*mu1;
		complex<double> Omega2 = 2.0 * (s+nu+mu2-mu1-1.0);
		complex<double> Omega3 = 1.0 - 2.0*s - 2.0*mu2;
		complex<double> Delta1 = s*(s+1.0) - l*(l+1.0) + 2.0*(s+mu2)*(mu1-nu)
									- nu + mu1 + mu2;
		complex<double> Delta2 = mu2 * (2.0*s+mu2);

		complex<double> alpha0 = Omega1 + 1.0, beta0 = Delta1;
		complex<double> xi1 = -beta0*xi0/alpha0;
		complex<double> xik = 0.0, xikm1 = xi1, xikm2 = xi0;

		complex<double> sum = xi0 + xi1*(r-1.0)/r;

		double r_coefficient = (r-1.0)/r;

		bool converged = false;
		for (int k = 2; k < kmax; ++k)
		{
			r_coefficient *= (r-1.0)/r;
			double km1 = k - 1;
			complex<double> alpha_km1 = 1.0 + (Omega1 + 1.0)/km1 + Omega1 / (km1*km1);
			complex<double> beta_km1 = -2.0 + (Omega2 + 2.0)/km1 + Delta1 / (km1*km1);
			complex<double> gamma_km1 = 1.0 - (Omega3 - 3.0)/km1 + (Delta2 - Omega3 + 2.0) / (km1*km1);
			
			xik = -( beta_km1 * xikm1 + gamma_km1 * xikm2 ) / alpha_km1;

			sum += xik * r_coefficient;
			if ( check_for_convergence and 2.0*abs(xik * r_coefficient) / abs(sum) < epsilon )
			{
				converged = true;
				break;	//converged; terminate loop
			}
			xikm2 = xikm1;
			xikm1 = xik;
		}

		sum *= prefactor;
		if (check_for_convergence and not converged)
			cerr << "Did not converge at r = " << r
					<< "!!  effective epsilon = " << 2.0*abs(xik * r_coefficient) / abs(sum)
					<< "; " << 2.0*abs(xik * r_coefficient) << "   " << abs(sum) << endl;
		cout << setw(16) << setprecision(12) << r << "   " << sum.real() << "   " << sum.imag() << endl;
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
