#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

int errfunc_integrand(double x, const double y[], double dydx[], void* params){
	dydx[0]=2.0/sqrt(M_PI)*exp(-x*x);
	return GSL_SUCCESS;
}

double myerf(double z){
	int arraydim = 1;
	gsl_odeiv2_system erfcsys;
	erfcsys.function = errfunc_integrand;
	erfcsys.dimension = arraydim;
	erfcsys.jacobian = NULL;
	erfcsys.params = NULL;


	double hstart, epsabs, epsrel;
	epsabs=epsrel=1e-6;
	hstart=copysign(0.1,z);

	gsl_odeiv2_driver* erfcdriver = gsl_odeiv2_driver_alloc_y_new(&erfcsys,
		gsl_odeiv2_step_rkf45, hstart, epsabs, epsrel);

	double y[1]={0}; //initial condition
	double a=0;

	gsl_odeiv2_driver_apply(erfcdriver,&a,z,y);

	gsl_odeiv2_driver_free(erfcdriver);

	return y[0];
	}

int main(int argc, char** argv){
assert(argc>=4);
	double a = atof(argv[1]); 
	double b = atof(argv[2]); 
	double dx = atof(argv[3]); 
for(double x=a;x<=b+dx;x+=dx)
	printf("%g %g\n",x,myerf(x));
return 0;
}
