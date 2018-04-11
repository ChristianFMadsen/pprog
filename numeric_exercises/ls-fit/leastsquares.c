#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "qr-ls.h"

double funs(int i, double x){
   switch(i){
   case 0: return log(x); break;
   case 1: return 1.0;   break;
   case 2: return x;     break;
   default: {fprintf(stderr,"funs: wrong i:%d",i); return NAN;}
   }
}


void leastsquares(double funs(int i, double x), int colDim, gsl_vector* x,  gsl_vector* y,  gsl_vector* errors,  gsl_vector* c,  gsl_matrix* S){

int rowDim = x->size;

gsl_matrix* I = gsl_matrix_alloc(colDim,colDim); gsl_matrix_set_identity(I);
gsl_matrix* A = gsl_matrix_alloc(rowDim,colDim);
gsl_vector* b = gsl_vector_alloc(rowDim);
gsl_matrix* R = gsl_matrix_alloc(colDim,colDim);
gsl_matrix* B = gsl_matrix_alloc(colDim,colDim);

for(int i=0; i<rowDim; i++){
	double x_i = gsl_vector_get(x,i);
	double y_i = gsl_vector_get(y,i);
	double errors_i = gsl_vector_get(errors,i);
	gsl_vector_set(b,i,y_i/errors_i); //b according to eq. (7)

	for(int k=0; k<colDim; k++){
		gsl_matrix_set(A, i, k, funs(k, x_i)/errors_i); } //A_ik according to eq. (7)
}

qr_gs_decomp(A,R); //A turns into Q
qr_gs_solve(A,R,b,c); //Solves QRc=b i.e. finds the vector c (eq 4) 
qr_gs_inverse(I,R,B); //Calculates inverse of A=QR and stores it in B s.t. AB=QRB=AA⁻¹=1
gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,B,B,0,S); //Calculates (A⁻¹)*(A⁻¹)^T (covariance matrix)


gsl_matrix_free(I); gsl_matrix_free(A); gsl_vector_free(b); gsl_matrix_free(R); gsl_matrix_free(B);
}