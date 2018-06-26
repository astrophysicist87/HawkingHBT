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
	compute_Rin();

	compute_Rout();

	compute_Rc();

	return 0;
}

