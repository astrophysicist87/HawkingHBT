#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

const int dim = 4;
const double l = 0.0;
const double omega = 0.4;

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
	 = gsl_odeiv_control_y_new (1e-6, 0.0);
	gsl_odeiv_evolve * e 
	 = gsl_odeiv_evolve_alloc (dim);

	double mu = 10;
	gsl_odeiv_system sys = {func, jac, dim, &mu};

	//double t = 1.0+1.e-8, t1 = 20.0;
	double t = 20.0, t1 = 1.0+1.e-8;
	double h = -1e-6;
	double y[dim] = { 1.0, 0.0, 0.0, 0.0 };

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
	double mu = *(double *)params;
	double c1r = 1.0 / (r*(r-1.0));
	double c1i = -2.0*omega*r*r / (r*(r-1.0));
	double c2 = -(l*(l+1.0) + 1.0/r) / (r*(r-1.0));

	//variable order: 0-Rr, 1-Ri, 2-Sr, 3-Si
	f[0] = y[2];
	f[1] = y[3];
	f[2] = -c1r * y[2] + c1i * y[3] - c2 * y[0];
	f[3] = -c1i * y[2] - c1r * y[3] - c2 * y[1];

	return GSL_SUCCESS;
}

int jac (double r, const double y[], double *dfdy, 
  double dfdt[], void *params)
{
	double mu = *(double *)params;
	double c1r = 1.0 / (r*(r-1.0));
	double c1i = -2.0*omega*r*r / (r*(r-1.0));
	double c2 = -(l*(l+1.0) + 1.0/r) / (r*(r-1.0));

	gsl_matrix_view dfdy_mat 
	 = gsl_matrix_view_array (dfdy, dim, dim);
	gsl_matrix * m = &dfdy_mat.matrix; 

	//initialize matrix to zero
	for (int i = 0; i < dim; ++i)
	for (int j = 0; j < dim; ++j)
		gsl_matrix_set (m, i, j, 0.0);

	gsl_matrix_set (m, 0, 2, 1.0);

	gsl_matrix_set (m, 1, 3, 1.0);

	gsl_matrix_set (m, 2, 0, -c2);
	gsl_matrix_set (m, 2, 2, -c1r);
	gsl_matrix_set (m, 2, 3, c1i);

	gsl_matrix_set (m, 3, 1, -c2);
	gsl_matrix_set (m, 3, 2, -c1i);
	gsl_matrix_set (m, 3, 3, -c1r);

	for (int i = 0; i < dim; ++i)
		dfdt[i] = 0.0;

	return GSL_SUCCESS;
}

