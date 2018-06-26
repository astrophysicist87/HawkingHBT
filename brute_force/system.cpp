#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <complex>

#include "RWeqnSolver.h"
#include "system.h"

using namespace std;


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



/////////////
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



