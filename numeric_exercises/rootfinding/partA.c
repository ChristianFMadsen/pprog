#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include "qr-rf.h"
#define eps 1e-6
#define dx 1e-8

//System of equations
int f1(gsl_vector* x1, gsl_vector* fx){
static int funCalls=0;
funCalls++;
double x=gsl_vector_get(x1,0), y=gsl_vector_get(x1,1);
int A=10000;
gsl_vector_set(fx,0,A*x*y-1.0);
gsl_vector_set(fx,1,exp(-x)+exp(-y)-1.0-1.0/A);
return funCalls;
}
//Rosenbrock
int f2(gsl_vector* x1, gsl_vector* fx){
static int funCalls=0;
funCalls++;
double x=gsl_vector_get(x1,0), y=gsl_vector_get(x1,1);
gsl_vector_set(fx,0,(-2)*(1-x)+2*100*(y-x*x)*(-2)*x);
gsl_vector_set(fx,1,2*100*(y-x*x));
return funCalls;
}
//Himmelblau
int f3(gsl_vector* x1, gsl_vector* fx){
static int funCalls=0;
funCalls++;
double x=gsl_vector_get(x1,0), y=gsl_vector_get(x1,1);
gsl_vector_set(fx,0,4*x*x*x+4*x*y-42*x+2*y*y-14);
gsl_vector_set(fx,1,2*x*x+4*x*y+4*y*y*y-26*y-22);
return funCalls;
}

int main() {
printf("Part A: Newton's method with back-tracking linesearch and numerical Jacobian\n");
int numberOfSteps, funCalls;
gsl_vector* x1=gsl_vector_alloc(2); gsl_vector* fx=gsl_vector_alloc(2);
printf("Finding the root of the system of equations:\n");
printf("A*x*y=1\n"); printf("exp(-x)+exp(-y)=1+1/A\n"); printf("where A=10000.\n");
printf("Start vector x:\n");
gsl_vector_set(x1,0,2); gsl_vector_set(x1,1,-2);
gsl_vector_fprintf(stdout, x1, "%g");
f1(x1,fx); 
printf("f(x):\n"); gsl_vector_fprintf(stdout, fx, "%g");
numberOfSteps=newton(f1,x1,dx,eps);
printf("Solution vector x:\n"); gsl_vector_fprintf(stdout, x1, "%g");
funCalls=f1(x1,fx);
printf("f(x):\n"); gsl_vector_fprintf(stdout, fx, "%g");
printf("Number of steps taken: %i\n", numberOfSteps);
printf("Number of function calls: %i\n", funCalls);

printf("\n"); printf("Finding the root of the gradient of the Rosenbrock function:\n"); 
printf("Start vector x:\n");
gsl_vector_set(x1,0,3); gsl_vector_set(x1,1,5);
gsl_vector_fprintf(stdout, x1, "%g");
f2(x1,fx); 
printf("f(x):\n"); gsl_vector_fprintf(stdout, fx, "%g");
numberOfSteps=newton(f2,x1,dx,eps);
printf("Solution vector x:\n"); gsl_vector_fprintf(stdout, x1, "%g");
funCalls=f2(x1,fx);
printf("f(x):\n"); gsl_vector_fprintf(stdout, fx, "%g");
printf("Number of steps taken: %i\n", numberOfSteps);
printf("Number of function calls: %i\n", funCalls);
printf("This yields f(1,1)=(1-1)²+100(1-1²)²=0 for the Rosenbrock function.\n");

printf("\n"); printf("Finding the root of the gradient of the Himmelblau function:\n"); 
printf("Start vector x:\n");
gsl_vector_set(x1,0,0); gsl_vector_set(x1,1,7);
gsl_vector_fprintf(stdout, x1, "%g");
f3(x1,fx); 
printf("f(x):\n"); gsl_vector_fprintf(stdout, fx, "%g");
numberOfSteps=newton(f3,x1,dx,eps);
printf("Solution vector x:\n"); gsl_vector_fprintf(stdout, x1, "%g");
funCalls=f3(x1,fx);
printf("f(x):\n"); gsl_vector_fprintf(stdout, fx, "%g");
printf("Number of steps taken: %i\n", numberOfSteps);
printf("Number of function calls: %i\n", funCalls);
printf("This yields f(3,2)=(3²+2-11)²+(3+2²-7)²=0 for the Himmelblau function.\n");
printf("\n");


gsl_vector_free(x1); gsl_vector_free(fx);
return 0;
}