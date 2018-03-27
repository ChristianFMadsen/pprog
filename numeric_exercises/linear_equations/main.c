#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "lin_eq.h"
#define RND ((double)rand()/RAND_MAX)
#define FORMAT "%8.3f"


void printMatrix(gsl_matrix* A){
	for(int i=0;i<A->size1;i++){
		for(int j=0;j<A->size2;j++) printf(FORMAT,gsl_matrix_get(A,i,j));
		printf("\n");}
}


int main(){
	size_t row_dim=8, col_dim=4;
	gsl_matrix* A = gsl_matrix_calloc(row_dim, col_dim);
	gsl_matrix* R = gsl_matrix_calloc(col_dim, col_dim);
	for (int i = 0; i < row_dim; i++)
	{
		for (int j = 0; j < col_dim; j++)
		{
			gsl_matrix_set(A,i,j,RND);
		}
	}
	printf("Modified GramSchmidt QR decomposition:\n");
	printf("Contents of randomly generated matrix A:\n");
	printMatrix(A);
	qr_gs_decomp(A,R);
	printf("Matrix Q:\n");
	printMatrix(A);
	printf("Matrix R:\n");
	printMatrix(R);
	printf("Note that R is upper triangular.\n");
	printf("Checking that Q^(TRANSPOSED)*Q=I:\n");
	gsl_matrix* QTQ = gsl_matrix_calloc(col_dim, col_dim);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, A, A, 0.0, QTQ);
	printMatrix(QTQ); printf("Success!\n");
	printf("Checking that QR=A:\n");
	gsl_matrix* QR = gsl_matrix_calloc(row_dim, col_dim);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, R, 0.0, QR);
	printMatrix(QR); printf("Success!\n"); printf("\n\n");

	printf("Random square matrix A:\n");
	gsl_matrix* A2 = gsl_matrix_calloc(col_dim, col_dim);
	gsl_matrix* Q2 = gsl_matrix_calloc(col_dim, col_dim);
	gsl_matrix* R2 = gsl_matrix_calloc(col_dim, col_dim);
	for (int i = 0; i < col_dim; i++)
	{
		for (int j = 0; j < col_dim; j++)
		{
			gsl_matrix_set(A2,i,j,RND);
		}
	}
	gsl_matrix_memcpy(Q2, A2);
	printMatrix(A2);
	printf("Random vector b:\n");
	gsl_vector* b = gsl_vector_calloc(col_dim);
	for (int i = 0; i < col_dim; i++)
	{
		gsl_vector_set(b, i, RND);
	}
	gsl_vector_fprintf(stdout, b, FORMAT);
	printf("Factorizing A=QR and solving QRx=b:\n");
	qr_gs_decomp(Q2,R2);
	gsl_vector* x = gsl_vector_calloc(col_dim);
	qr_gs_solve(Q2, R2, b, x);
	printf("Computed x:\n");
	gsl_vector_fprintf(stdout, x, FORMAT);
	printf("Computed Q:\n");
	printMatrix(Q2);
	printf("Computed R:\n");
	printMatrix(R2);
	printf("Checking that Ax=b:\n");
	gsl_blas_dgemv(CblasNoTrans, 1.0, A2, x, 0.0, b);
	gsl_vector_fprintf(stdout, b, FORMAT); printf("Success!\n");

	gsl_matrix* A_inv = gsl_matrix_calloc(col_dim, col_dim);
	qr_gs_inverse(Q2, R2, A_inv);
	printf("Computed inverse of A ie. A^(-1):\n"); 
	printMatrix(A_inv);
	printf("Checking that A*A^(-1)=I:\n");
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A2,A_inv,0.0,R2);
	printMatrix(R2); printf("Success!\n");
	printf("Checking that A^(-1)*A=I:\n"); 
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A_inv,A2,0.0,R2);
	printMatrix(R2); printf("Success!\n");

	gsl_matrix_free(A); gsl_matrix_free(R); gsl_matrix_free(A2); gsl_matrix_free(Q2); gsl_matrix_free(R2); gsl_vector_free(b); gsl_vector_free(x); gsl_matrix_free(A_inv); gsl_matrix_free(QTQ); gsl_matrix_free(QR);
	return 0;
}