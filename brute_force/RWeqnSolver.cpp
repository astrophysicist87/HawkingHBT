#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <complex>

#include "RWeqnSolver.h"

using namespace std;

int func (double t, const double y[], double f[],
   void *params);
int jac (double t, const double y[], double *dfdy, 
  double dfdt[], void *params);

int main (void)
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

	//double t = 1.0+1.e-8, t1 = 20.0;
	double t = rInf, t1 = r0;
	double h = -1e-6;
	double rsInf = rInf + 2.0*M*log( rInf/(2.0*M) - 1.0 );
	//do purely ingoing mode as test (will not be same as EF solution)
	complex<double> Rin_at_rInf = exp(-i * omega * rsInf) / rInf;
	complex<double> Rin_prime_at_rInf
						= - exp(-i * omega * rsInf)
							* ( 1.0/rInf + i*omega*rInf/(rInf-2.0*M) )
							/ rInf;
	double y[dim] = { Rin_at_rInf.real(),
						Rin_at_rInf.imag(),
						Rin_prime_at_rInf.real(),
						Rin_prime_at_rInf.imag() };

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

	return 0;
}
int func (double r, const double y[], double f[],
	void *params)
{
	double local_Veff = Veff(l, M, r);
	double a = 1.0/r + 1.0/(r-2.0*M);
	double b = ( r*r*r*r*(omega*omega - local_Veff) + 2.0*M*r - 4.0*M*M )
				/ ( r*r*(r-2.0*M)*(r-2.0*M) );

	//variable order: 0-Rr, 1-Ri, 2-Sr, 3-Si
	f[0] = y[2];
	f[1] = y[3];
	f[2] = -a * y[2] - b * y[0];
	f[3] = -a * y[3] - b * y[1];

	return GSL_SUCCESS;
}

int jac (double r, const double y[], double *dfdy, 
	double dfdt[], void *params)
{
	double local_Veff = Veff(l, M, r);
	double a = 1.0/r + 1.0/(r-2.0*M);
	double b = ( r*r*r*r*(omega*omega - local_Veff) + 2.0*M*r - 4.0*M*M )
				/ ( r*r*(r-2.0*M)*(r-2.0*M) );

	gsl_matrix_view dfdy_mat 
	 = gsl_matrix_view_array (dfdy, dim, dim);
	gsl_matrix * m = &dfdy_mat.matrix; 

	//initialize matrix to zero
	for (int i = 0; i < dim; ++i)
	for (int j = 0; j < dim; ++j)
		gsl_matrix_set (m, i, j, 0.0);
	gsl_matrix_set (m, 0, 2, 1.0);

	gsl_matrix_set (m, 1, 3, 1.0);

	gsl_matrix_set (m, 2, 0, -b);
	gsl_matrix_set (m, 2, 2, -a);

	gsl_matrix_set (m, 3, 1, -b);
	gsl_matrix_set (m, 3, 3, -a);

	for (int i = 0; i < dim; ++i)
		dfdt[i] = 0.0;

	return GSL_SUCCESS;
}



