#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "adapInt.h"

int fCalls; //Global variable
double int5GSL(double x, void* params){
	fCalls++;
	return x/(exp(x)-1); 
}

int main(){
printf("Testing recursive adaptive integrator:\n"); 
double a=0, b=1, err, acc=1e-6, eps=1e-6;
printf("Integrating sqrt(x) on the interval [%g,%g] with acc=%g and eps=%g\n", a,b,acc,eps);
double Q = my_integrator(int1, a, b, acc, eps, &err,&fCalls);
printf("integral_result=%lg, estimated error=%lg\n\n",Q,err);
printf("Integrating 1/sqrt(x) on the interval [%g,%g] with acc=%g and eps=%g\n", a,b,acc,eps);
Q=my_integrator(int2,a,b,acc,eps,&err,&fCalls);
printf("integral_result=%lg, estimated error=%lg\n\n",Q,err);
printf("Integrating ln(x)/sqrt(x) on the interval [%g,%g] with acc=%g and eps=%g\n", a,b,acc,eps);
Q=my_integrator(int3,a,b,acc,eps,&err,&fCalls);
printf("integral_result=%lg, estimated error=%lg\n\n",Q,err);
fCalls=0;
acc=1e-11; eps=1e-10;
printf("Integrating 4*sqrt(1-(1-x)^2) on the interval [%g,%g] with acc=%g and eps=%g\n", a,b,acc,eps);
Q=my_integrator(int4,a,b,acc,eps,&err,&fCalls);
printf("integral_result=%.50lg, estimated error=%lg, function calls:%i\n",Q,err,fCalls);
double diff=M_PI-Q;
printf("Difference between math.h M_PI and the integral result: %.50lg \n\n", diff);
acc=1e-10; eps=1e-10; fCalls=0;
printf("Integrating x/(exp(x)-1) on the interval [%g,Inf) with acc=%g and eps=%g\n",a,acc,eps);
printf("Result: pi^2/6 = 1.64493...\n");
Q=aToInf(int5,a,acc,eps,&err,&fCalls);
printf("integral_result=%g, estimated error=%g, function calls: %i\n\n", Q,err,fCalls);	
printf("Using GSL routines:\n");
gsl_integration_workspace* w=gsl_integration_workspace_alloc(1000); fCalls=0;
gsl_function F; F.function=&int5GSL; F.params=NULL;
gsl_integration_qagiu(&F,a,acc,eps,1000,w,&Q,&err);
printf("integral_result=%g, error=%g, function calls: %i\n\n", Q,err,fCalls);
a=-1; b=1;
printf("Integrating sin(x)/sqrt(x+1) on the interval [%g,%g] with acc=%g and eps=%g\n", a,b,acc,eps);
printf("Note that the integrand is divergent at x=-1.\n");
printf("Result=-0.8266196541509...\n");
printf("First using the integrator from part A:\n");
fCalls=0;
Q=my_integrator(int6,-1,1,acc,eps,&err,&fCalls);
printf("integral_result=%.13g, estimated error=%g, function calls: %i\n\n", Q, err, fCalls);
printf("Now using the Clenshaw-Curtis transformation:\n");
fCalls=0;
Q=clenCurtis(int6,acc,eps,&err,&fCalls);
printf("integral_result=%.13g, estimated error=%g, function calls: %i\n\n", Q, err, fCalls);
printf("Note that the Clenshaw-Curtis transformation reduces the amount of function calls and produces a result closer to the actual result compared to the integrator from part A.");

gsl_integration_workspace_free(w);
return 0;
}