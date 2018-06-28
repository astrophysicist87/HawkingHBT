#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <complex>

#include "RWeqnSolver.h"
#include "Rout.h"
#include "system.h"

using namespace std;

void compute_Rout()
{
	const gsl_odeiv_step_type * T 
	 = gsl_odeiv_step_rk8pd;
	//const gsl_odeiv_step_type * T 
	// = gsl_odeiv_step_rk4imp;

	gsl_odeiv_step * s 
	 = gsl_odeiv_step_alloc (T, dim);
	gsl_odeiv_control * c 
	 = gsl_odeiv_control_y_new (relative_precision, 0.0);
	gsl_odeiv_evolve * e 
	 = gsl_odeiv_evolve_alloc (dim);

	double M_loc = M;
	gsl_odeiv_system sys = {func, jac, dim, &M_loc};

	double t = rInf, t1 = r0;
	double h = -1e-6;
	double rsInf = ( M < 1.e-6 ) ? rInf : rInf + 2.0*M*log( rInf/(2.0*M) - 1.0 );
	//do purely ingoing mode as test (will not be same as EF solution)
	complex<double> Rout_at_rInf = exp(i * omega * rsInf) / rInf;
	complex<double> Rout_prime_at_rInf
						= - exp(i * omega * rsInf)
							* ( 1.0/rInf - i*omega*rInf/(rInf-2.0*M) )
							/ rInf;
	double y[dim] = { Rout_at_rInf.real(),
						Rout_at_rInf.imag(),
						Rout_prime_at_rInf.real(),
						Rout_prime_at_rInf.imag() };

	while (t > t1)
	{
		int status = gsl_odeiv_evolve_apply (e, c, s,
					                        &sys, 
					                        &t, t1,
					                        &h, y);

		if (status != GSL_SUCCESS)
			break;

		//printf ("ROUT: %.5e %.5e %.5e %.5e %.5e\n", t, y[0], y[1], y[2], y[3]);
		//printf ("ROUT: %.5e %.15e %.15e\n", t, y[0], y[1]);
	}

	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);

	return;
}


