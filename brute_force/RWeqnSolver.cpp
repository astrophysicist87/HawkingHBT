#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <complex>

#include "RWeqnSolver.h"
#include "Rin.h"
#include "Rout.h"
#include "Rc.h"

using namespace std;

int main (void)
{
	//compute_Rin();

	//compute_Rout();

	double rsInf = rInf + 2.0*M*log( rInf/(2.0*M) - 1.0 );
	complex<double> Rc_at_rInf(0,0), Rc_prime_at_rInf(0,0);
	complex<double> Rin_at_rInf = exp(-i * omega * rsInf) / rInf;
	complex<double> Rin_prime_at_rInf
						= - exp(-i * omega * rsInf)
							* ( 1.0/rInf + i*omega*rInf/(rInf-2.0*M) )
							/ rInf;
	complex<double> Rout_at_rInf = exp(i * omega * rsInf) / rInf;
	complex<double> Rout_prime_at_rInf
						= - exp(i * omega * rsInf)
							* ( 1.0/rInf - i*omega*rInf/(rInf-2.0*M) )
							/ rInf;

	//get the causal solution
	compute_Rc( Rc_at_rInf, Rc_prime_at_rInf );

	//get the coefficients for in and out modes
	complex<double> Ain(0,0), Aout(0,0);
	solve_2D_simult_eqns( Rin_at_rInf, Rin_prime_at_rInf,
							Rout_at_rInf, Rout_prime_at_rInf,
							Rc_at_rInf, Rc_prime_at_rInf,
							Ain, Aout );

	printf ("Ain = %.15e, Aout = %.15e\n", Ain, Aout);

	return 0;
}

