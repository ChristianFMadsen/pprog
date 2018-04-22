#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <math.h>
#include "minimization.h"

int main(){
	int n=2; //Number of variables (x and y)
	double eps=1e-5, dx=1e-5;
	printf("Part B: Quasi-Newton method with Broyden's update:\n");
	printf("Finding the minimum of the Rosenbrock valley function: \n");
	printf("Start vector x:\n");
	gsl_vector* xstart=gsl_vector_alloc(n);
	gsl_vector_set(xstart,0,3); gsl_vector_set(xstart,1,5);
	gsl_vector_fprintf(stdout,xstart,"%g");
	int numberOfSteps;
	numberOfSteps=qNewton(fx_ros,xstart,eps,dx);
	printf("Minimum found at:\n");
	gsl_vector_fprintf(stdout,xstart,"%g");
	printf("At found minimum: f(x)=%g\n", fx_ros(xstart));
	printf("In %i steps.\n", numberOfSteps);

	printf("Finding the minimum of the Himmelblau function:\n");
	printf("Start vector x:\n");
	gsl_vector_set(xstart,0,-1); gsl_vector_set(xstart,1,7);
	gsl_vector_fprintf(stdout,xstart,"%g");
	numberOfSteps=qNewton(fx_him,xstart,eps,dx);
	printf("Minimum found at:\n");
	gsl_vector_fprintf(stdout,xstart,"%g");
	printf("At found minimum: f(x)=%g\n", fx_him(xstart));
	printf("In %i steps.\n", numberOfSteps);	

	printf("Compared to the root finding exercise where the Newton's method with numerical Jacobian was used:\n"); 
	printf("Rosenbrock: 393 steps taken, Himmelblau: 8 steps taken\n");
	//Part B.iv)
	double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
	double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
	double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
	int N = sizeof(t)/sizeof(t[0]);
	
	FILE* datastream=fopen("data.txt","w");
	fprintf(datastream, "# t[i] \t y[i] \t e[i]\n");
	for (int i = 0; i < N; i++){
		fprintf(datastream, "%g \t %g \t %g \n", t[i], y[i], e[i]); }

	gsl_vector* xstart_decay=gsl_vector_alloc(3);
	gsl_vector_set(xstart_decay, 0, -5); gsl_vector_set(xstart_decay, 1, 5); gsl_vector_set(xstart_decay, 2, 10);
	qNewton(decay_fun,xstart_decay,eps,dx);
	double A=gsl_vector_get(xstart_decay,0); double T=gsl_vector_get(xstart_decay,1); double B=gsl_vector_get(xstart_decay,2);
	
	fprintf(datastream, "\n\n");
	fprintf(datastream, "# t \t A*exp(-t/T)+B\n");
	for (double t = 0; t < 10; t=t+0.1)
	{
		fprintf(datastream, "%g \t %g\n", t, A*exp(-t/T)+B);
	}

	gsl_vector_free(xstart); gsl_vector_free(xstart_decay);
	fclose(datastream);
	return 0;
}