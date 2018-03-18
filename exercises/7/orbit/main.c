#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

int eq1(double x, const double y[], double dydx[], void* params){
	dydx[0]=y[0]*(1-y[0]);
	return GSL_SUCCESS;
}

int orbitalfunc(double phi, const double u[], double dudphi[], void* params){
	double epsilon=*(double*)params;
	dudphi[0]=u[1];
	dudphi[1]=1-u[0]+epsilon*u[0]*u[0];
	return GSL_SUCCESS;
}

int main(){
	//Part 1 of "Problems orbit"
	int arraydim=1;
	gsl_odeiv2_system sys;
	sys.function = eq1;
	sys.jacobian = NULL;
	sys.dimension = arraydim;
	sys.params = NULL;

	gsl_odeiv2_driver* d=gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 1e-6);
	double xini=0.0, xfin=3.0; //xinitial and xfinal
	double y[1]={0.5}; //initial value
	double dx=0.05;


	for(double i=xini; i<xfin; i=i+dx){
	int status = gsl_odeiv2_driver_apply(d,&xini,i,y);
	if(status != GSL_SUCCESS){
		printf("error, return value=%i\n", status);
		break;
	}


	printf("%g \t %g\n",i,y[0]);
	}

	printf("\n"); printf("\n");

	//Part 2 of "problems orbit":
	int arraydimOrb=2;
	double epsilon;
	gsl_odeiv2_system sys2;
	sys2.function = orbitalfunc;
	sys2.jacobian = NULL;
	sys2.dimension = arraydimOrb;
	sys2.params = (void*) &epsilon;



	gsl_odeiv2_driver* dOrb=gsl_odeiv2_driver_alloc_y_new(&sys2, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 1e-6); //Last three args: initial stepsize, epsabs, epsrel
	
	double u[2];

	double initvals[3]={0.0, -0.5, -0.5};
	double epsilonvals[3]={0.0, 0.0, 0.01};


	for (int k = 0; k <=2; k++){
		double xiniOrb=0.0, xfinOrb=196*M_PI;
		epsilon=epsilonvals[k];
		u[0]=1;
		u[1]=initvals[k];
			for(double phi=xiniOrb; phi<xfinOrb; phi=phi+dx){
				int status = gsl_odeiv2_driver_apply(dOrb, &xiniOrb, phi,u);
					if(status != GSL_SUCCESS){
					printf("error, return value=%i\n", status);
					break;
					}

			printf("%g \t %g\n",phi,u[0]);
			}

		printf("\n"); printf("\n");
	}

	gsl_odeiv2_driver_free(d);
	gsl_odeiv2_driver_free(dOrb);
	return 0;
}