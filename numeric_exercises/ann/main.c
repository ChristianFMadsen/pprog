#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include "ann.h"
#define RND (double)rand()/RAND_MAX

double activation_function(double x){
	return x*exp(-x*x); }

double fitFun(double x){
	return cos(x)*sin(x)*cos(2*x)*sin(2*x)*exp(-x*x*x); }

double fitFun2d(double x1, double x2){
	return x1*x1-x2*x2; } 

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

	printf("\n\n"); //Part B. 
	n=10; a=-5; b=5; nx=30;  
	ann2d* network2d=ann_alloc2d(n,activation_function);
	gsl_matrix* xys=gsl_matrix_alloc(nx*nx,2); 
 	gsl_vector* flist=gsl_vector_alloc(nx*nx);
	for(int i=0; i<nx; i++){ //Print tabulated x,y and f(x,y) values. (index 2)
		for(int j=0; j<nx; j++){
		double x1=a+(b-a)*i/(nx-1); //x
		double x2=a+(b-a)*j/(nx-1); //y
		double f=fitFun2d(x1,x2);   //z
		gsl_matrix_set(xys,i*nx+j,0,x1);
		gsl_matrix_set(xys,i*nx+j,1,x2);
		gsl_vector_set(flist,i*nx+j,f);	
		printf("%g %g %g\n", x1, x2, f);
		}
	}

	for(int i=0;i<network2d->n;i++){ //Provide network with data
		gsl_vector_set(network2d->data,0*network2d->n+i,a+(b-a)*i/(network2d->n-1));
		gsl_vector_set(network2d->data,1*network2d->n+i,1);
		gsl_vector_set(network2d->data,2*network2d->n+i,1);
		gsl_vector_set(network2d->data,3*network2d->n+i,a+(b-a)*i/(network2d->n-1));
		gsl_vector_set(network2d->data,4*network2d->n+i,1);
		gsl_vector_set(network2d->data,5*network2d->n+i,1); }

	ann_train2d(network2d,xys,flist);

	printf("\n\n"); 
	for(int i=0; i<nx; i++){ //Print out fit data. (index 3)
		for(int j=0; j<nx; j++){
		double x1=a+(b-a)*i/(nx-1); 					//x
		double x2=a+(b-a)*j/(nx-1); 					//y
		double f=ann2d_feed_forward(network2d,x1,x2);	//z
		printf("%g %g %g\n", x1, x2, f);
		}
	}

ann_free(network); ann2d_free(network2d); gsl_vector_free(vx); gsl_vector_free(vf); gsl_matrix_free(xys); gsl_vector_free(flist);
return 0;
}