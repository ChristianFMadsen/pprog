#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include "qr-rf.h"
/* void f(gsl_vector* x, gsl_vector* fx) contains a guess in x and then calculates the function value f(x) */
int newton(int f(gsl_vector* x1, gsl_vector* fx), gsl_vector* x, double dx, double eps){
int n=x->size;
int numberOfSteps=0;
gsl_matrix* J = gsl_matrix_alloc(n,n); //Jacobian matrix
gsl_matrix* R = gsl_matrix_alloc(n,n); //R matrix from QR decomp of J.
gsl_vector* fx = gsl_vector_alloc(n); //vector that holds f(x) i.e. the functional value at x
gsl_vector* z  = gsl_vector_alloc(n); //used to update x with
gsl_vector* fz = gsl_vector_alloc(n); //used to update f(x) with	
gsl_vector* df = gsl_vector_alloc(n); //holds f(x+dx) when building Jacobian
gsl_vector* Dx = gsl_vector_alloc(n); // $\Delta \vec{x}$ 

while(1){
//Build Jacobian matrix
f(x,fx);
for(int k=0; k<n; k++){
	gsl_vector_set(x,k,gsl_vector_get(x,k)+dx); //update x
	f(x,df); //calculate f(x+dx) and store in df
	gsl_vector_sub(df,fx); //calculate df-f(x)=f(x+dx)-f(x) and store it in df
	for(int i=0; i<n; i++){ 
		gsl_matrix_set(J,i,k,gsl_vector_get(df,i)/dx); } //Update J_ik'th entry of jacobian matrix with finite difference derivative eq. (7)
	gsl_vector_set(x,k,gsl_vector_get(x,k)-dx); //subtract dx to preserve v
}
qr_gs_decomp(J,R); //QR decomp: J -> Q
qr_gs_solve(J,R,fx,Dx); //Solve QR(Dx)=fx i.e. finds the vector Dx (eq. 5 up to a sign)
gsl_vector_scale(Dx, -1); //Fixes the sign

double lambda=1;
while(1){
	gsl_vector_memcpy(z,x); //copy vec x into vec z
	gsl_vector_add(z,Dx); //x+Delta_x
	f(z,fz); //calculate f(x+Delta_x) and store in fz

	/* The following is the bit of pseudocode under eq. 8: */
	if(gsl_blas_dnrm2(fz)<(1-lambda/2)*gsl_blas_dnrm2(fx) || lambda<0.015625){
		break; } //dnrm2 calculates euclidean norm, 1/64=0.015625.
	lambda=lambda*0.5;
	gsl_vector_scale(Dx,0.5);
	}
/*Update x and f(x) even if the condition eq. 8 couldn't be met. 
If it is not met the step is taken away from this difficult place
 with the minimum lambda value s.t. it can try again */
gsl_vector_memcpy(x,z); //Update x s.t. x is now x+Delta_x
gsl_vector_memcpy(fx,fz); //Update f(x) s.t. f(x) is now f(x+Delta_x)
numberOfSteps++;
/*Terminate search if the norm of Delta_x becomes smaller than dx 
or the norm of f(x) becomes smaller than epsilon 
i.e. it has converged to the necessary accuracy */
if( gsl_blas_dnrm2(Dx)<dx || gsl_blas_dnrm2(fx)<eps ){ 
	break;
	}
}

gsl_matrix_free(J); gsl_matrix_free(R); gsl_vector_free(fx); gsl_vector_free(fz); gsl_vector_free(z); gsl_vector_free(df); gsl_vector_free(Dx);
return numberOfSteps;
}


int newton_with_jacobian(int f(const gsl_vector* x, gsl_vector* fx, gsl_matrix* J), gsl_vector* x, double eps){
int n=x->size;
int numberOfSteps=0;
gsl_matrix* J = gsl_matrix_alloc(n,n);
gsl_matrix* fJ = gsl_matrix_alloc(n,n); 
gsl_matrix* R = gsl_matrix_alloc(n,n); 
gsl_vector* fx = gsl_vector_alloc(n); 
gsl_vector* z  = gsl_vector_alloc(n); 
gsl_vector* fz = gsl_vector_alloc(n); 
gsl_vector* Dx = gsl_vector_alloc(n); 

while(1){
f(x,fx,J);
qr_gs_decomp(J,R);
qr_gs_solve(J,R,fx,Dx); 
gsl_vector_scale(Dx, -1); 

double lambda=1;
while(1){
	gsl_vector_memcpy(z,x); 
	gsl_vector_add(z,Dx); 
	f(z,fz,fJ);
	if(gsl_blas_dnrm2(fz)<(1-lambda/2)*gsl_blas_dnrm2(fx) || lambda<0.015625){
		break; } //dnrm2 calculates euclidean norm, 1/64=0.015625.
	lambda=lambda*0.5;
	gsl_vector_scale(Dx,0.5);
	}
gsl_vector_memcpy(x,z);
gsl_vector_memcpy(fx,fz); 
numberOfSteps++;

if(gsl_blas_dnrm2(fx)<eps){ 
	break; }
}

gsl_matrix_free(J); gsl_matrix_free(fJ); gsl_matrix_free(R); gsl_vector_free(fx); gsl_vector_free(z); gsl_vector_free(fz); gsl_vector_free(Dx); 
return numberOfSteps;
}

int newton_with_jacobian_refined(int f(const gsl_vector* x, gsl_vector* fx, gsl_matrix* J), gsl_vector* x, double eps){
int n=x->size;
int numberOfSteps=0;
gsl_matrix* J = gsl_matrix_alloc(n,n); 
gsl_matrix* fJ = gsl_matrix_alloc(n,n);
gsl_matrix* R = gsl_matrix_alloc(n,n); 
gsl_vector* fx = gsl_vector_alloc(n); 
gsl_vector* z  = gsl_vector_alloc(n); 
gsl_vector* fz = gsl_vector_alloc(n); 
gsl_vector* Dx = gsl_vector_alloc(n); 

do{
	f(x,fx,J);
	qr_gs_decomp(J,R);
	qr_gs_solve(J,R,fx,Dx); 
	gsl_vector_scale(Dx, -1); 
	double g0 = 0.5*pow(gsl_blas_dnrm2(fx),2);
	double gprime0 = -pow(gsl_blas_dnrm2(fx),2);
	double glambda = 0;
	double lambda=1;
	do{
		double c = (glambda-g0-gprime0*lambda)/(lambda*lambda);
		glambda = g0+gprime0*lambda+c*lambda*lambda;
		gsl_vector_memcpy(z,x);
		gsl_vector_scale(Dx,lambda);
		gsl_vector_add(z,Dx);
		f(z,fz,fJ);
		lambda=gprime0/(2*c);
	}while(gsl_blas_dnrm2(fz) > (1-lambda/2)*gsl_blas_dnrm2(fx) && lambda > 1.0/64);
	gsl_vector_memcpy(x,z);
	gsl_vector_memcpy(fx,fz);
	numberOfSteps++;
}while(gsl_blas_dnrm2(fx)>eps);

gsl_matrix_free(J); gsl_matrix_free(fJ); gsl_matrix_free(R); gsl_vector_free(fx); gsl_vector_free(z); gsl_vector_free(fz); gsl_vector_free(Dx); 
return numberOfSteps;
}