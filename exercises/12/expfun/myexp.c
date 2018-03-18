#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_vector.h>

int exp_integrand(double x, const double y[], double dydx[], void* params)
{
	dydx[0]=y[0];
	return GSL_SUCCESS;
}


double myexp(double x)
{
if(x<0){x=(-1)*x; return 1/myexp(x);}
if(x>=1) return myexp(x/2)*myexp(x/2);
assert(x>=0 && x<1);

int arraydim = 1;
gsl_odeiv2_system expsys;
expsys.function = exp_integrand;
expsys.dimension = arraydim;
expsys.jacobian = NULL;
expsys.params = NULL;


double hstart, epsabs, epsrel;
epsabs=0;
epsrel=1e-10;
hstart=copysign(1e-6,x);

gsl_odeiv2_driver* expdriver = gsl_odeiv2_driver_alloc_y_new(&expsys,
	gsl_odeiv2_step_rkf45, hstart, epsabs, epsrel);


double y[1]={1}; //initial condition
double startPoint=0;

int status = gsl_odeiv2_driver_apply(expdriver,&startPoint,x,y);

if (status!=GSL_SUCCESS) printf("error, return value=%i\n", status); 

gsl_odeiv2_driver_free(expdriver);

return y[0];
}
