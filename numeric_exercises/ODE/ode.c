#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <assert.h>
#include "ode.h"

//y''=-y
void sinODE(double x, gsl_vector* y, gsl_vector* dydx){
    double y1 = gsl_vector_get(y,0), y2 = gsl_vector_get(y,1);
    gsl_vector_set(dydx,0,y2); gsl_vector_set(dydx,1,-y1);
}

double intFun(double x){  
    return exp(x)*x*x*x*sin(x*x);
}


void rkstep12(double x, double h, gsl_vector* yx, void f(double x, gsl_vector* y, gsl_vector* dydx), gsl_vector* yxh, gsl_vector* err){
int n=yx->size;
gsl_vector* k0=gsl_vector_alloc(n);
gsl_vector* yt=gsl_vector_alloc(n);
gsl_vector* k12=gsl_vector_alloc(n);
f(x,yx,k0); // calculate k0 (eq. 14)
for(int i=0;i<n;i++) { // y value to be used in calculating k12 (eq. 14)
	gsl_vector_set(yt,i,gsl_vector_get(yx,i)+gsl_vector_get(k0,i)*h/2); }

f(x+h/2,yt,k12); // calculate k1/2 (eq. 14)
for(int i=0;i<n;i++) { // calculate y_(i+1) (eq. 12)
	gsl_vector_set(yxh,i,gsl_vector_get(yx,i)+gsl_vector_get(k12,i)*h); }

// error estimate (eq. 23)
for(int i=0;i<n;i++){
	gsl_vector_set(err,i,(gsl_vector_get(k0,i)-gsl_vector_get(k12,i))*h/2); } 

gsl_vector_free(k0); gsl_vector_free(yt); gsl_vector_free(k12);
}


int ode_driver(gsl_vector* xlist, gsl_matrix* ylist, double b, double h, double acc, double eps, 
	int maxSteps, int dim, void f(double x, gsl_vector* y,gsl_vector* dydx)){

int n=dim;
int numberOfSteps=0; //int i; 
double x, errNorm, yNorm, tol, a=gsl_vector_get(xlist,0);
assert(a!=b);

gsl_vector* yx=gsl_vector_alloc(n);
gsl_vector* yxh=gsl_vector_alloc(n);
gsl_vector* err=gsl_vector_alloc(n);
while(gsl_vector_get(xlist,numberOfSteps)<b){ //Keep stepping until x=b   
x=gsl_vector_get(xlist,numberOfSteps); gsl_matrix_get_row(yx,ylist,numberOfSteps); //current x value and y values. 
if(x + h > b){
h = b-x; } //modify step size if the step size makes x overshoot the end point 

rkstep12(x,h,yx,f,yxh,err); //Use stepper rkstep12
 
errNorm=gsl_blas_dnrm2(err); //calculate local error (eq. 41)
yNorm=gsl_blas_dnrm2(yxh);

tol=(yNorm*eps+acc)*sqrt(h/(b-a)); // calculate local tolerance (eq. 41)
if(errNorm<tol){ // if err<tol accept the step
numberOfSteps++; 

	if(numberOfSteps>=maxSteps){
		printf("ERROR: Maximum number of steps reached without converging.\n");
		return -numberOfSteps; }

gsl_vector_set(xlist,numberOfSteps,x + h); //Save x-step
for(int i=0;i<n;i++){ //Save y-step
gsl_matrix_set(ylist,numberOfSteps,i,gsl_vector_get(yxh,i)); }
}


if(err>0){ //update step size (eq. 40)
	h *= pow(tol/errNorm,0.25)*0.95; }
else{ //If err=0 try a bigger step size
	h *= 2;
} 
	 	
}/*end while*/
return numberOfSteps;

gsl_vector_free(yx); gsl_vector_free(yxh); gsl_vector_free(err);
}

double integrator(gsl_vector* xlist, gsl_matrix* ylist, double b, double h, double acc, double eps, int maxSteps, int dim, double f(double x)){
    
    void integSystem(double x, gsl_vector* y, gsl_vector* dydx){
        double fx=f(x);
        gsl_vector_set(dydx,0,fx); }

    int numberOfSteps=ode_driver(xlist,ylist,b,h,acc,eps,maxSteps,dim,integSystem);
    //printf("number of steps taken is: %i\n",k);
    return gsl_matrix_get(ylist,numberOfSteps,0);
}



