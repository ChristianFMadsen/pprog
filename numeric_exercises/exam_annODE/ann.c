#include <math.h>
#include "ann.h"

ann* ann_alloc(int n, double (*g)(double), double (*dg)(double)){
	ann* network = malloc(sizeof(ann));
	network->n=n;
	network->g=g;
	network->dg=dg;
	network->data=gsl_vector_alloc(3*n);
	return network;
}


void ann_free(ann* network){
	gsl_vector_free(network->data);
	free(network);
}


void ann_feed_forward(ann* network, double x, double* F, double* dF){ 
	int size=network->n;
	double Fres=0, dFres=0;
	for(int i=0; i<size; i++){
		double a=gsl_vector_get(network->data,0*size+i); 
		double b=gsl_vector_get(network->data,1*size+i); 
		double w=gsl_vector_get(network->data,2*size+i); 
		Fres += network->g((x-a)/b)*w; 
		dFres += network->dg((x-a)/b)*w/b;
	}
	*F = Fres;
	*dF = dFres;
}


void ann_train(ann* network, gsl_vector* vx, double (*func)(double, double), double x0, double y0){
	double dx=1e-6, eps=1e-5;
	int N=vx->size;
	double Fpxk, dFpxk, Fpx0, dFpx0;
	double deviation(gsl_vector* p){ //Deviation function
		gsl_vector_memcpy(network->data,p);
		double sum=0;
		for(int k=0; k<N; k++){
			double xk = gsl_vector_get(vx,k);
			ann_feed_forward(network, xk, &Fpxk, &dFpxk);
			double fk = func(xk, Fpxk);
			sum += pow(fabs((dFpxk-fk)),2); //Calculates the sum in the deviation function
		}
		ann_feed_forward(network,x0,&Fpx0, &dFpx0);
		sum+=N*pow(fabs(Fpx0-y0),2); //Adds the last term

		return sum/N;
	}
	gsl_vector* p=gsl_vector_alloc(network->data->size);
	gsl_vector_memcpy(p,network->data);
	qNewton(deviation, p, eps, dx); //Minimize deviation function
	gsl_vector_memcpy(network->data,p);
	gsl_vector_free(p);
}