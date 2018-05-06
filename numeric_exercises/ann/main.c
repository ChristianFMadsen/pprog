#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include "ann.h"
double activation_function(double x){
	return x*exp(-x*x); }

double fitFun(double x){
	return cos(x)*sin(x)*cos(2*x)*sin(2*x)*exp(-x*x*x); }

int main(){
	int n=5;
	ann* network=ann_alloc(n,activation_function);
	double a=-1,b=1; //fit from a to b
	int nx=20;
	gsl_vector* vx=gsl_vector_alloc(nx); //tabulated x values
	gsl_vector* vf=gsl_vector_alloc(nx); //tabulated f(x) values
	for(int i=0;i<nx;i++){
		double x=a+(b-a)*i/(nx-1); //Evenly spaced x's from a to b
		double f=fitFun(x); 
		gsl_vector_set(vx,i,x); 
		gsl_vector_set(vf,i,f); }

	for(int i=0;i<network->n;i++){ //Provide network with data
		gsl_vector_set(network->data,0*network->n+i,a+(b-a)*i/(network->n-1));
		gsl_vector_set(network->data,1*network->n+i,1);
		gsl_vector_set(network->data,2*network->n+i,1); }

	ann_train(network,vx,vf); //Train network

	for(int i=0;i<vx->size;i++){ //Print out the tabulated x and f(x) values.
		double x=gsl_vector_get(vx,i);
		double f=gsl_vector_get(vf,i);
		printf("%g %g\n",x,f);
	}
	printf("\n\n"); //fit data will be at index 1

	double dz=1.0/64;
	for(double z=a;z<=b;z+=dz){ //Print out data for the fitted curve.
		double y=ann_feed_forward(network,z);
		printf("%g %g\n",z,y);
	}
	
ann_free(network);
return 0;
}
