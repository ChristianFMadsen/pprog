#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include "ann.h"

double activation_function(double x){
	return x*exp(-x*x);
}

double activation_function_derivative(double x){
	return (1-2*x*x)*exp(-x*x);
}

double logistic_function(double x, double y){
	return y*(1-y);
}

double logistic_function_solution(double x){
	return 1/(1+exp(-x));
}

double gauss_function(double x, double y){
	return -x*y;
}

double gauss_function_solution(double x){
	return exp(-x*x/2.0);
}

int main(){
	FILE* examprob = fopen("ExamProblemAndDetails.txt", "w");
	fprintf(examprob, "This is a solution of the exam problem 'Artificial neural network (ANN) for solving ODE'.\nMade by Christian F. Madsen, student number: 201506198\n");
	int N=20;
	int n=7; //amount of hidden neurons
	double a=-5.0, b=5.0, x0=0.0, y0Log=0.5, y0Gauss=1.0; //start and end points as well as initial conditions

	gsl_vector* vx=gsl_vector_alloc(N); //vector with x-values spanning [a,b]
	for(int k=1; k<=N; k++){
		double xk=a+(b-a)*(k-1)/(N-1);
		gsl_vector_set(vx,k-1,xk);
	}

	ann* network_logistic = ann_alloc(n,&activation_function,&activation_function_derivative);  //network to solve the logistic function ODE
	for(int i=0; i<n; i++){
		gsl_vector_set(network_logistic->data,0*network_logistic->n+i,a+(b-a)*i/(network_logistic->n-1));
		gsl_vector_set(network_logistic->data,1*network_logistic->n+i,1);
		gsl_vector_set(network_logistic->data,2*network_logistic->n+i,1);
	}

	ann* network_gauss = ann_alloc(n,&activation_function,&activation_function_derivative); //network to solve the gaussian function ODE
	for(int i=0; i<n; i++){
		gsl_vector_set(network_gauss->data,0*network_gauss->n+i,a+(b-a)*i/(network_gauss->n-1));
		gsl_vector_set(network_gauss->data,1*network_gauss->n+i,1);
		gsl_vector_set(network_gauss->data,2*network_gauss->n+i,1);
	}

	ann_train(network_logistic, vx, &logistic_function, x0, y0Log); //train networks
	ann_train(network_gauss, vx, &gauss_function, x0, y0Gauss);

	FILE* logisticData=fopen("logisticData.txt","w"); 
	FILE* gaussData=fopen("gaussData.txt","w");

	double dx=10.0/100, logistic_ann, gauss_ann, dF;
	for(double xi=a; xi<b; xi+=dx){ //print ANN data as well as the analytical solution
		ann_feed_forward(network_logistic,xi,&logistic_ann,&dF);
		ann_feed_forward(network_gauss,xi,&gauss_ann,&dF);
		fprintf(logisticData, "%g \t %g \t %g\n", xi, logistic_ann, logistic_function_solution(xi));
		fprintf(gaussData, "%g \t %g \t %g\n", xi, gauss_ann, gauss_function_solution(xi));
	}
		
	fclose(examprob); fclose(logisticData); fclose(gaussData); gsl_vector_free(vx); ann_free(network_gauss); ann_free(network_logistic);
}