#ifndef SYSTEM_H
#define SYSTEM_H

#include <complex>

using namespace std;

int func (double t, const double y[], double f[],
   void *params);
int jac (double t, const double y[], double *dfdy, 
  double dfdt[], void *params);

#endif
