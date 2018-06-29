#include <math.h>
#include <stdlib.h>
#include "montecarlo.h"
#define RND ((double)rand()/RAND_MAX)

double int1(double* x){
	return tan(*x); 
}

double int2(double* x){
	return exp(pow(*x,2));
}

double int3(double* x){
	double x1, y, z;
	x1=x[0]; y=x[1]; z=x[2];
	return pow(1-cos(x1)*cos(y)*cos(z),-1)*pow(M_PI,-3);
}

void rndx(int dim, double* a, double* b, double* x){
for(int i=0; i<dim; i++){
	x[i]=a[i]+RND*(b[i]-a[i]); } //generates a random x in the integration interval since 0<=RND<=1 
}

void plainMC(int dim, double* a, double* b, double f(double* x), int N, double* result, double* error){ 
double V=1; 
for(int i=0; i<dim; i++){
	V*=b[i]-a[i]; //Calculate V^(dim) i.e. the auxiliary rectangular area/volume/...	
} 

double sum=0, sum2=0, fx, x[dim]; //x[dim] holds the dim x values in each direction.

for(int i=0;i<N;i++){ 
	rndx(dim,a,b,x); //generate random x 
	fx=f(x); //Calculate integrand at this x
	sum+=fx; 
	sum2+=fx*fx; }

double avg=sum/N, var=sum2/N-avg*avg; // Calculate average and sigma^2 (variance) according to eq. 4 and 2
*result=avg*V; //Calculate result according to eq. 2
*error=sqrt(var/N)*V; //Calculate error according to eq. 3
}



double int2D(double f(double* x, int* fCalls), int* fCalls, double a, double b, double c(double), double d(double), double acc, double eps, double* error){
	double f1(double x1, int* fCalls){
		double f2(double x2, int* fCalls){
			double xs[]={x1,x2};
			return f(xs, fCalls);
		}
		return my_integrator(f2, c(x1), d(x1), acc, eps, error, fCalls);
	}
	return my_integrator(f1, a, b, acc, eps, error, fCalls);
}