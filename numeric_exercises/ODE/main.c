#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <math.h>
#include "ode.h"

int main() {
FILE* datastream = fopen("sin.txt","w");

printf("Solving the ODE representation of the sine function:\n"); printf("See the result in sinPlot.pdf\n");
double b=5*M_PI, h=0.1, acc=1e-4, eps=0;
int numberOfSteps, dim=2, maxSteps=1e5, vecDim=10000;
gsl_vector* xlist = gsl_vector_alloc(vecDim);
gsl_vector_set(xlist,0,0);
gsl_matrix* ylist = gsl_matrix_alloc(vecDim,dim);
gsl_matrix_set(ylist,0,0,0); gsl_matrix_set(ylist,0,1,1);

numberOfSteps=ode_driver(xlist,ylist,b,h,acc,eps,maxSteps,dim,sinODE);
printf("Steps taken: %i\n", numberOfSteps);
fprintf(datastream, "#x \t y(x) \t y'(x) \n");
for(int i=0; i<numberOfSteps; i++){
 	fprintf(datastream, "%g %g %g\n", gsl_vector_get(xlist,i), gsl_matrix_get(ylist,i,0), gsl_matrix_get(ylist,i,1));
}

printf("Solving the integral: exp(x)*x^3*sin(x^2) from x=0 to x=exp(1)\n");
dim=1; h=1e-3; b=exp(1); eps=1e-3;
gsl_vector* integxlist=gsl_vector_alloc(vecDim); gsl_matrix* integylist=gsl_matrix_alloc(vecDim,dim);
gsl_vector_set(integxlist,0,0);
gsl_matrix_set(integylist,0,0,0);
double integResult=integrator(integxlist,integylist,b,h,acc,eps,maxSteps,dim,intFun);
printf("Analytical result: -9.00154\n");
printf("ODE integrator result: %g \n", integResult);

fclose(datastream); gsl_matrix_free(ylist); gsl_vector_free(xlist);
return 0;
}