#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include "minimization.h"

int main(){
	int n=2; //Number of variables (x and y)
	double eps=1e-5;
	printf("Part A: Newton's minimization method with explicit derivatives:\n");
	printf("Finding the minimum of the Rosenbrock valley function: \n");
	printf("Start vector x:\n");
	gsl_vector* xstart=gsl_vector_alloc(n);
	gsl_vector_set(xstart,0,3); gsl_vector_set(xstart,1,5);
	gsl_vector_fprintf(stdout,xstart,"%g");
	int numberOfSteps;
	numberOfSteps=newton(f_ros,xstart,eps);
	printf("Minimum found at:\n");
	gsl_vector_fprintf(stdout,xstart,"%g");
	printf("At found minimum: f(x)=%g\n", fx_ros(xstart));
	printf("In %i steps.\n", numberOfSteps);

	printf("Finding the minimum of the Himmelblau function:\n");
	printf("Start vector x:\n");
	gsl_vector_set(xstart,0,-1); gsl_vector_set(xstart,1,7);
	gsl_vector_fprintf(stdout,xstart,"%g");
	numberOfSteps=newton(f_him,xstart,eps);
	printf("Minimum found at:\n");
	gsl_vector_fprintf(stdout,xstart,"%g");
	printf("At found minimum: f(x)=%g\n", fx_him(xstart));
	printf("In %i steps.\n", numberOfSteps);	

	printf("\n\n");
	gsl_vector_free(xstart);
	return 0;
}