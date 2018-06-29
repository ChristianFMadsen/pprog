#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include "qr-rf.h"
#define eps 1e-6

//System of equations
int f1(const gsl_vector* x1, gsl_vector* fx, gsl_matrix* J){
static int funCalls=0;
funCalls++;
double x=gsl_vector_get(x1,0), y=gsl_vector_get(x1,1);
int A=10000;
gsl_vector_set(fx,0,A*x*y-1.0);
gsl_vector_set(fx,1,exp(-x)+exp(-y)-1.0-1.0/A);
gsl_matrix_set(J,0,0,A*y); gsl_matrix_set(J,0,1,A*x);
gsl_matrix_set(J,1,0,-exp(-x)); gsl_matrix_set(J,1,1,-exp(-y));
return funCalls;
}

//Rosenbrock
int f2(const gsl_vector* x1, gsl_vector* fx, gsl_matrix* J){
static int funCalls=0;
funCalls++;
double x=gsl_vector_get(x1,0), y=gsl_vector_get(x1,1);
gsl_vector_set(fx,0,(-2)*(1-x)+2*100*(y-x*x)*(-2)*x);
gsl_vector_set(fx,1,2*100*(y-x*x));
gsl_matrix_set(J,0,0,1200*x*x-400*y+2); gsl_matrix_set(J,0,1,-400*x);
gsl_matrix_set(J,1,0,-400*x); gsl_matrix_set(J,1,1,200);
return funCalls;
}
//Himmelblau
int f3(const gsl_vector* x1, gsl_vector* fx, gsl_matrix* J){
static int funCalls=0;
funCalls++;
double x=gsl_vector_get(x1,0), y=gsl_vector_get(x1,1);
gsl_vector_set(fx,0,4*x*x*x+4*x*y-42*x+2*y*y-14);
gsl_vector_set(fx,1,2*x*x+4*x*y+4*y*y*y-26*y-22);
gsl_matrix_set(J,0,0,12*x*x+4*y-42); gsl_matrix_set(J,0,1,4*(x+y));
gsl_matrix_set(J,1,0,4*(x+y)); gsl_matrix_set(J,1,1,12*y*y+4*x-26);
return funCalls;
}

int main() {
printf("\n\nPart C: Newton's method with analytic Jacobian and refined linesearch.\n");
int numberOfSteps, funCalls;
gsl_vector* x1=gsl_vector_alloc(2); 
gsl_vector* fx=gsl_vector_alloc(2); 
gsl_matrix* J=gsl_matrix_alloc(2,2);
printf("Finding the root of the system of equations:\n");
printf("A*x*y=1\n"); printf("exp(-x)+exp(-y)=1+1/A\n"); printf("where A=10000.\n");
printf("Start vector x:\n");
gsl_vector_set(x1,0,2); gsl_vector_set(x1,1,-2);
gsl_vector_fprintf(stdout, x1, "%g");
f1(x1,fx,J); 
printf("f(x):\n"); gsl_vector_fprintf(stdout, fx, "%g");
numberOfSteps=newton_with_jacobian_refined(f1,x1,eps);
printf("Solution vector x:\n"); gsl_vector_fprintf(stdout, x1, "%g");
funCalls=f1(x1,fx,J);
printf("f(x):\n"); gsl_vector_fprintf(stdout, fx, "%g");
printf("Number of steps taken: %i\n", numberOfSteps);
printf("Number of function calls: %i\n", funCalls);

printf("\n"); printf("Finding the root of the gradient of the Rosenbrock function:\n"); 
printf("Start vector x:\n");
gsl_vector_set(x1,0,3); gsl_vector_set(x1,1,5);
gsl_vector_fprintf(stdout, x1, "%g");
f2(x1,fx,J); 
printf("f(x):\n"); gsl_vector_fprintf(stdout, fx, "%g");
numberOfSteps=newton_with_jacobian_refined(f2,x1,eps);
printf("Solution vector x:\n"); gsl_vector_fprintf(stdout, x1, "%g");
funCalls=f2(x1,fx,J);
printf("f(x):\n"); gsl_vector_fprintf(stdout, fx, "%g");
printf("Number of steps taken: %i\n", numberOfSteps);
printf("Number of function calls: %i\n", funCalls);
printf("This yields f(1,1)=(1-1)²+100(1-1²)²=0 for the Rosenbrock function.\n");

printf("\n"); printf("Finding the root of the gradient of the Himmelblau function:\n"); 
printf("Start vector x:\n");
gsl_vector_set(x1,0,0); gsl_vector_set(x1,1,7);
gsl_vector_fprintf(stdout, x1, "%g");
f3(x1,fx,J); 
printf("f(x):\n"); gsl_vector_fprintf(stdout, fx, "%g");
numberOfSteps=newton_with_jacobian_refined(f3,x1,eps);
printf("Solution vector x:\n"); gsl_vector_fprintf(stdout, x1, "%g");
funCalls=f3(x1,fx,J);
printf("f(x):\n"); gsl_vector_fprintf(stdout, fx, "%g");
printf("Number of steps taken: %i\n", numberOfSteps);
printf("Number of function calls: %i\n", funCalls);
printf("This yields f(3,2)=(3²+2-11)²+(3+2²-7)²=0 for the Himmelblau function.\n");
printf("\n");


printf("Comparisons:\n");
printf("The number of steps taken does not change between the method with a numerical Jacobian and the method with an analytic Jacobian.\n");
printf("The amount of function calls is reduced when using the method with an analytic Jacobian however.\n");
printf("In all cases the number of steps and function calls is reduced when using the method with an analytic Jacobian combined with a refined linesearch. \n");
printf("The method with an analytic Jacobian combined with the refined linesearch even takes the same number of steps as the GSL routine.\n");

gsl_vector_free(x1); gsl_vector_free(fx); gsl_matrix_free(J); 
return 0;
}