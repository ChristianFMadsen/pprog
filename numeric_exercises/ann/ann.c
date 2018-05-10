#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <stdio.h>
#include <math.h>
#include "ann.h"

ann* ann_alloc(int n, double (*f)(double x)){
	ann* network = malloc(sizeof(ann));
	network->n=n;
	network->f=f;
	network->data=gsl_vector_alloc(3*n);
	return network;
}


void ann_free(ann* network){
	gsl_vector_free(network->data);
	free(network);
}


double ann_feed_forward(ann* network, double x){ //calculates output of hidden layer neurons
	double sum=0, size=network->n;
	for(int i=0; i<size; i++){
		double a=gsl_vector_get(network->data,0*size+i); //first n entries is the a_i's, next n entries is the b_i's, and last n entries is the w_i's.
		double b=gsl_vector_get(network->data,1*size+i);
		double w=gsl_vector_get(network->data,2*size+i);
		sum=sum+network->f((x+a)/b)*w;
	}
	return sum;
}


void ann_train(ann* network, gsl_vector* vx, gsl_vector* vf){
	double dx=1e-9, eps=1e-5;
	double deviation(gsl_vector* p){ 
		gsl_vector_memcpy(network->data,p);
		double sum=0;
		for(int i=0; i<vx->size; i++){
			double x=gsl_vector_get(vx,i); //tabulated x-values of function to be interpolated
			double f=gsl_vector_get(vf,i); //tabulated function values of function to be interpolated
			double y=ann_feed_forward(network,x); 
			sum=sum+pow(y-f,2);
		}
		return sum/vx->size;
	}
	gsl_vector* p=gsl_vector_alloc(network->data->size);
	gsl_vector_memcpy(p,network->data);
	qNewton(deviation, p, eps, dx);
	gsl_vector_memcpy(network->data,p);
	gsl_vector_free(p);

}





ann2d* ann_alloc2d(int n, double (*f)(double x)){
	ann2d* network = malloc(sizeof(ann2d));
	network->n=n;
	network->f=f;
	network->data=gsl_vector_alloc(6*n); 
	return network;
}


void ann2d_free(ann2d* network){
	gsl_vector_free(network->data);
	free(network);
}

void ann_train2d(ann2d* network, gsl_matrix* xys, gsl_vector* vf){
	double dx=1e-9, eps=1e-1;
	double deviation(gsl_vector* p){ 
		gsl_vector_memcpy(network->data,p);
		double sum=0;
		for(int i=0; i<vf->size; i++){
			double x1=gsl_matrix_get(xys,i,0);
			double x2=gsl_matrix_get(xys,i,1);
			double f=gsl_vector_get(vf,i); //tabulated function values of function to be interpolated
			double fAnn=ann2d_feed_forward(network,x1,x2); 
			sum=sum+pow(f-fAnn,2);
		}
		return sum/vf->size;
	}
	gsl_vector* p=gsl_vector_alloc(network->data->size);
	gsl_vector_memcpy(p,network->data);
	qNewton(deviation, p, eps, dx);
	gsl_vector_memcpy(network->data,p);
	gsl_vector_free(p);

}

double ann2d_feed_forward(ann2d* network, double x1, double x2){
	double sum=0, arg1, arg2; 
	int size=network->n;
	for(int i=0; i<size; i++){
		double a1=gsl_vector_get(network->data,0*size+i);
		double b1=gsl_vector_get(network->data,1*size+i);
		double w1=gsl_vector_get(network->data,2*size+i);
		double a2=gsl_vector_get(network->data,3*size+i);
		double b2=gsl_vector_get(network->data,4*size+i);
		double w2=gsl_vector_get(network->data,5*size+i);
		arg1=(x1+a1)/b1; 
		arg2=(x2+a2)/b2;
		sum=sum+network->f(arg1)*w1+network->f(arg2)*w2;
	}
	return sum;
}