#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include "ann.h"

void gradient(double f(gsl_vector* x), gsl_vector* x, gsl_vector* df, double dx){
int n=x->size; 
double dfi, fx, fxdx, xi;
fx=f(x);
for(int i=0; i<n; i++){
		xi = gsl_vector_get(x,i);
		gsl_vector_set(x,i,xi+dx);
		fxdx=f(x);
		gsl_vector_set(x,i,xi);
		dfi=(fxdx-fx)/dx;
		gsl_vector_set(df,i,dfi);
	}
}

double f_ros(gsl_vector* x1, gsl_vector* df, gsl_matrix* H){
	double x=gsl_vector_get(x1,0); double y=gsl_vector_get(x1,1);
	gsl_vector_set(df,0,2*x-400*x*(y-x*x)-2);
	gsl_vector_set(df,1,200*y-200*x*x);
	gsl_matrix_set(H,0,0,1200*x*x-400*y+2);
	gsl_matrix_set(H,0,1,-400*x);
	gsl_matrix_set(H,1,0,-400*x);
	gsl_matrix_set(H,1,1,200);
	return (1-x)*(1-x)+100*(y-x*x)*(y-x*x);
}

double fx_ros(gsl_vector* x1){
	double x=gsl_vector_get(x1,0); double y=gsl_vector_get(x1,1);
	return (1-x)*(1-x)+100*(y-x*x)*(y-x*x);
}

double f_him(gsl_vector* x1, gsl_vector* df, gsl_matrix* H){
	double x=gsl_vector_get(x1,0); double y=gsl_vector_get(x1,1);
	gsl_vector_set(df,0,4*x*x*x+4*x*y-42*x+2*y*y-14);
	gsl_vector_set(df,1,2*x*x+4*x*y+4*y*y*y-26*y-22);
	gsl_matrix_set(H,0,0,12*x*x+4*y-42);
	gsl_matrix_set(H,0,1,4*(x+y));
	gsl_matrix_set(H,1,0,4*(x+y));
	gsl_matrix_set(H,1,1,12*y*y+4*x-26);
	return (x*x+y-11)*(x*x+y-11)+(x+y*y-7)*(x+y*y-7);
}

double fx_him(gsl_vector* x1){
	double x=gsl_vector_get(x1,0); double y=gsl_vector_get(x1,1);
	return (x*x+y-11)*(x*x+y-11)+(x+y*y-7)*(x+y*y-7);
}


double decay_fun(gsl_vector* x){
	double  A = gsl_vector_get(x,0);
	double  T = gsl_vector_get(x,1);
	double  B = gsl_vector_get(x,2);
	double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
	double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
	double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
	int n = sizeof(t)/sizeof(t[0]);
	double sum=0;
	for(int i=0;i<n;i++)
	{
	 sum = sum + pow((A*exp(-t[i]/T)+B - y[i])/e[i],2);
	}
	return sum;
}


//xstart: starting point, becomes the latest approximation to the root on exit 
//double f(...): f: objective function, df: gradient, H: Hessian matrix H
//eps: accuracy goal, on exit |gradient|<eps
int newton(double f(gsl_vector* x, gsl_vector* df, gsl_matrix* H), gsl_vector* xstart, double eps){
int n = xstart->size;
int numberOfSteps=0;

gsl_vector* x=gsl_vector_calloc(n);
gsl_vector* z=gsl_vector_calloc(n);
gsl_vector* lambdaDeltax=gsl_vector_calloc(n); 
gsl_vector* deltax=gsl_vector_calloc(n);
gsl_vector* df=gsl_vector_calloc(n);
gsl_matrix* H=gsl_matrix_calloc(n,n);
gsl_matrix* R=gsl_matrix_calloc(n,n);
gsl_matrix* invH=gsl_matrix_calloc(n,n);

gsl_vector_set(x,0,gsl_vector_get(xstart,0));
gsl_vector_set(x,1,gsl_vector_get(xstart,1));

while(1){
double fx=f(x, df, H); //f(x)
qr_gs_decomp(H,R); //H->Q
qr_gs_inverse(H,R,invH);
gsl_blas_dgemv(CblasNoTrans,-1.0,invH,df,0,deltax);


double lambda=1;
double alpha=0.01;
double dotProd; //Stores scalar product between delta x transposed and the gradient of f(x) (eq. 8)
gsl_blas_ddot(deltax, df, &dotProd);

while(1){
gsl_vector_memcpy(z,x);
gsl_vector_memcpy(lambdaDeltax,deltax);
gsl_vector_scale(lambdaDeltax,lambda); //lambda*Delta_x
gsl_vector_add(z,lambdaDeltax); //x+lambda*Delta_x
double fxdelx=f(z,df,H); //f(x+lambda*Delta_x)
double armijo=fx+alpha*lambda*dotProd;

if(fxdelx < armijo || lambda<0.015625){ //  1/64=0.015625
	break;
}
lambda=lambda*0.5;
}
gsl_vector_memcpy(x,z); //Updates x
numberOfSteps++;
f(x, df, H); //recalculate gradient i.e. df
if(gsl_blas_dnrm2(df)<eps){
	break; }
}

gsl_vector_memcpy(xstart,x);
gsl_vector_free(x); gsl_vector_free(z); gsl_vector_free(df); gsl_vector_free(deltax); gsl_matrix_free(H); gsl_vector_free(lambdaDeltax); gsl_matrix_free(R); gsl_matrix_free(invH);
return numberOfSteps;
}


int qNewton(double f(gsl_vector* x), gsl_vector* xstart, double eps, double dx){
int n=xstart->size;

gsl_vector* x=gsl_vector_alloc(n);
gsl_vector* dfdx=gsl_vector_alloc(n);
gsl_matrix* invH=gsl_matrix_alloc(n,n);
gsl_vector* Dx=gsl_vector_alloc(n);
gsl_vector* s=gsl_vector_alloc(n);
gsl_vector* z=gsl_vector_alloc(n);
gsl_vector* dfdx_z=gsl_vector_alloc(n);
gsl_matrix* y=gsl_matrix_alloc(n,1);
gsl_matrix* invHy=gsl_matrix_alloc(n,1);
gsl_matrix* sMatrix=gsl_matrix_alloc(n,1);
gsl_matrix* sMinusinvHy=gsl_matrix_alloc(n,1);
gsl_matrix* sTransinvH=gsl_matrix_alloc(1,n);
gsl_matrix* numerator=gsl_matrix_alloc(n,n);
gsl_vector* yVec=gsl_vector_alloc(n);
gsl_vector* invHs=gsl_vector_alloc(n);
int numberOfSteps=0;

for(int j=0; j<n; j++){ //x=xstart
gsl_vector_set(x,j,gsl_vector_get(xstart,j));	
}
gsl_matrix_set_identity(invH); //H⁻¹=I

gradient(f,x,dfdx,dx); //dfdx=grad f
double fx=f(x), alpha=0.1, fz, sTransdfdx, sTrans_s, denominator;

do{
gsl_blas_dgemv(CblasNoTrans,-1.0,invH,dfdx,0,Dx); //Dx=\delta x
gsl_vector_memcpy(s,Dx);
gsl_vector_scale(s,2.0);


while(1){
	gsl_vector_scale(s,0.5);
	gsl_vector_memcpy(z,s);
	gsl_vector_add(z,x);
	fz=f(z);
	gsl_blas_ddot(s,dfdx,&sTransdfdx);
	if(fz < fx+alpha*sTransdfdx){
		break;
	}

	gsl_blas_ddot(s,s,&sTrans_s);
	if(sTrans_s<dx){
		gsl_matrix_set_identity(invH); //reset inverse Hessian
		break;
	}
}

gradient(f,z,dfdx_z,dx);
for(int j=0; j<n; j++){ //y as given below eq. 11
	gsl_matrix_set(y,j,0,gsl_vector_get(dfdx_z,j)-gsl_vector_get(dfdx,j));
}

for(int j=0; j<n; j++){ //build s matrix
	gsl_matrix_set(sMatrix,j,0,gsl_vector_get(s,j));
}

gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,invH,y,0.0,invHy);
gsl_matrix_memcpy(sMinusinvHy,sMatrix);
gsl_matrix_sub(sMinusinvHy,invHy);
gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,sMatrix,invH,0.0,sTransinvH);
gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,sMinusinvHy,sTransinvH,0.0,numerator);

for(int j=0; j<n; j++){ //build y in vector format
	gsl_vector_set(yVec,j,gsl_matrix_get(y,j,0));
}
gsl_blas_dgemv(CblasNoTrans,1.0,invH,s,0.0,invHs);
gsl_blas_ddot(yVec,invHs,&denominator);
gsl_matrix_scale(numerator,1.0/denominator);
gsl_matrix_add(invH,numerator);


fx=fz;
gsl_vector_memcpy(x,z); numberOfSteps++;
gsl_vector_memcpy(dfdx,dfdx_z);
}while(gsl_blas_dnrm2(Dx)>dx && gsl_blas_dnrm2(dfdx)>eps);

gsl_vector_memcpy(xstart,x);

gsl_vector_free(x); gsl_vector_free(dfdx); gsl_vector_free(Dx); gsl_vector_free(s); gsl_vector_free(z); gsl_vector_free(dfdx_z); gsl_vector_free(yVec);
gsl_vector_free(invHs); gsl_matrix_free(invH); gsl_matrix_free(y); gsl_matrix_free(invHy); gsl_matrix_free(sMatrix); gsl_matrix_free(sMinusinvHy);
gsl_matrix_free(sTransinvH); gsl_matrix_free(numerator);
return numberOfSteps;
}