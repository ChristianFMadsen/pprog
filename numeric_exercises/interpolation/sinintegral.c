#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

double f (double x, void * params) {
  double f = sin(5*x);
  return f;
}

double sinintegral(double a, double b)
{

  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

  double result, epsabs, epsrel, abserr;

  epsabs=1e-8;
  epsrel=0;


  gsl_function F;
  F.function = &f;
  F.params = NULL;

  gsl_integration_qag(&F, a, b, epsabs, epsrel, 1000, 6, w, &result, 
    &abserr);


  gsl_integration_workspace_free (w);

  return result;
}