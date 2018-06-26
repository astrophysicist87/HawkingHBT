#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <complex>

#include "RWeqnSolver.h"
#include "Rin.h"
#include "Rout.h"

using namespace std;

int main (void)
{
	compute_Rin();

	compute_Rout();

	return 0;
}

