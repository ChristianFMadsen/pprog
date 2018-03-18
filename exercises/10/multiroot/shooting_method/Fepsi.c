#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>
#include<stdio.h>
#include<math.h>

int ode_H_s_wave(double r, const double y[], double dy[], void* params)
{
	double eps = *(double*) params;
	dy[0]=y[1];
	dy[1]=-2*(1/r+eps)*y[0];
return GSL_SUCCESS;
}

double Fepsi(double eps, double r)
{
	int probDim=2;
	const double rmin = 1e-3;
	if(r<rmin)
	{
		return r-r*r;
	}

	gsl_odeiv2_system sys;
	sys.function=ode_H_s_wave;
	sys.jacobian = NULL;
	sys.dimension = probDim;
	sys.params = (void*)&eps;

	double hstart, abserr, relerr;
	hstart=1e-3;
	abserr=1e-7;
	relerr=1e-7;

	gsl_odeiv2_driver* D = gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,hstart,abserr,relerr);

	double r1=rmin, y[]={r1-r1*r1, 1-2*r1};
	int status = gsl_odeiv2_driver_apply(D, &r1, r, y);
	if (status != GSL_SUCCESS) 
	{
		fprintf(stderr, "Fepsi: odeiv2 error: %i\n", status);
	}

gsl_odeiv2_driver_free(D);
return y[0];
}