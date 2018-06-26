#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <complex>

#include "RWeqnSolver.h"
#include "Rc.h"
#include "system.h"

using namespace std;

void compute_Rc()
{
	const gsl_odeiv_step_type * T 
	 = gsl_odeiv_step_rk8pd;

	gsl_odeiv_step * s 
	 = gsl_odeiv_step_alloc (T, dim);
	gsl_odeiv_control * c 
	 = gsl_odeiv_control_y_new (1e-10, 0.0);
	gsl_odeiv_evolve * e 
	 = gsl_odeiv_evolve_alloc (dim);

	double M_loc = M;
	gsl_odeiv_system sys = {func, jac, dim, &M_loc};

	double t = r0, t1 = rInf;
	double h = 1e-6;
	double rs0 = r0 + 2.0*M*log( r0/(2.0*M) - 1.0 );
	//do purely ingoing mode as test (will not be same as EF solution)
	complex<double> Rc_at_r0 = exp(i * omega * rsInf) / r0;
	complex<double> Rc_prime_at_r0
						= - exp(i * omega * rs0)
							* ( 1.0/r0 - i*omega*r0/(r0-2.0*M) )
							/ r0;
	double y[dim] = { Rc_at_r0.real(),
						Rc_at_r0.imag(),
						Rc_prime_at_r0.real(),
						Rc_prime_at_r0.imag() };

	while (t >= t1)
	{
		int status = gsl_odeiv_evolve_apply (e, c, s,
					                        &sys, 
					                        &t, t1,
					                        &h, y);

		if (status != GSL_SUCCESS)
			break;

		printf ("%.5e %.5e %.5e %.5e %.5e\n", t, y[0], y[1], y[2], y[3]);
	}

	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);

	return;
}


