#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <omp.h>
#include "adapInt.h"

int fCalls; //Global variable
double int5GSL(double x, void* params){
	fCalls++;
	return x/(exp(x)-1); 
}

int main(){
printf("Using multiprocessing to evaluate integrals using my adaptive integrator.\n");
printf("Thread 1 will handle part A of the exercise while thread 2 will handle part B and C of the exercise:\n");

#pragma omp parallel sections
{
	#pragma omp section
	{
	printf("----[1] First thread started on part A:----\n");
	printf("[1] Part A: Testing recursive adaptive integrator:\n"); 
	double a=0, b=1, err, acc=1e-6, eps=1e-6;
	printf("[1] Integrating sqrt(x) on the interval [%g,%g] with acc=%g and eps=%g\n", a,b,acc,eps);
	double Q = my_integrator(int1, a, b, acc, eps, &err,&fCalls);
	printf("[1] integral_result=%lg, estimated error=%lg\n\n",Q,err);
	printf("[1] Integrating 1/sqrt(x) on the interval [%g,%g] with acc=%g and eps=%g\n", a,b,acc,eps);
	Q=my_integrator(int2,a,b,acc,eps,&err,&fCalls);
	printf("[1] integral_result=%lg, estimated error=%lg\n\n",Q,err);
	printf("[1] Integrating ln(x)/sqrt(x) on the interval [%g,%g] with acc=%g and eps=%g\n", a,b,acc,eps);
	Q=my_integrator(int3,a,b,acc,eps,&err,&fCalls);
	printf("[1] integral_result=%lg, estimated error=%lg\n\n",Q,err);
	fCalls=0;
	acc=1e-11; eps=1e-10;
	printf("[1] Integrating 4*sqrt(1-(1-x)^2) on the interval [%g,%g] with acc=%g and eps=%g\n", a,b,acc,eps);
	Q=my_integrator(int4,a,b,acc,eps,&err,&fCalls);
	printf("[1] integral_result=%.50lg, estimated error=%lg, function calls:%i\n",Q,err,fCalls);
	double diff=M_PI-Q;
	printf("[1] Difference between math.h M_PI and the integral result: %.50lg \n\n", diff);
	printf("----[1] First thread ended on part A.----\n");
	}

	#pragma omp section
	{
	printf("----[2] Second thread started on part B:----\n");
	double acc=1e-10, eps=1e-10, a=0, b, err; fCalls=0;
	printf("[2] Part B: Integrating x/(exp(x)-1) on the interval [%g,Inf) with acc=%g and eps=%g\n",a,acc,eps);
	printf("[2] Result: pi^2/6 = 1.64493...\n");
	double Q=aToInf(int5,a,acc,eps,&err,&fCalls);
	printf("[2] integral_result=%g, estimated error=%g, function calls: %i\n\n", Q,err,fCalls);	
	printf("[2] Using GSL routines:\n");
	gsl_integration_workspace* w=gsl_integration_workspace_alloc(1000); fCalls=0;
	gsl_function F; F.function=&int5GSL; F.params=NULL;
	gsl_integration_qagiu(&F,a,acc,eps,1000,w,&Q,&err);
	printf("[2] integral_result=%g, error=%g, function calls: %i\n\n", Q,err,fCalls);
	printf("[2] Part C: Clenshaw-Curtis variable transformation:\n");
	a=-1; b=1;
	printf("[2] Integrating sin(x)/sqrt(x+1) on the interval [%g,%g] with acc=%g and eps=%g\n", a,b,acc,eps);
	printf("[2] Note that the integrand is divergent at x=-1.\n");
	printf("[2] Result=-0.8266196541509...\n");
	printf("[2] First using the integrator from part A:\n");
	fCalls=0;
	Q=my_integrator(int6,-1,1,acc,eps,&err,&fCalls);
	printf("[2] integral_result=%.13g, estimated error=%g, function calls: %i\n\n", Q, err, fCalls);
	printf("[2] Now using the Clenshaw-Curtis transformation:\n");
	fCalls=0;
	a=-1; b=1;
	Q=clenCurtis(int6,a,b,acc,eps,&err,&fCalls);
	printf("[2] integral_result=%.13g, estimated error=%g, function calls: %i\n\n", Q, err, fCalls);
	printf("[2] Note that the Clenshaw-Curtis transformation reduces the amount of function calls and produces a result closer to the actual result compared to the integrator from part A.\n");
	gsl_integration_workspace_free(w);
	printf("----[2] Second thread ended on part B.----\n");
	}
}

return 0;
}