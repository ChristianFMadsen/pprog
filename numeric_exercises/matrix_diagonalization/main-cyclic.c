#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "md.h"
#define RND ((double)rand()/RAND_MAX)
#define FORMAT "%8.3f"


int main(int argc, char** argv)
{
	assert(argc>1);
	int dim = atoi(argv[1]);
	printf("Matrix dimension=%i\n", dim);

	gsl_matrix* A = gsl_matrix_calloc(dim,dim);
	gsl_matrix* B = gsl_matrix_calloc(dim,dim);

	for (int i = 0; i < dim; i++){
		for (int j = i; j < dim; j++){ 
			double rnum = RND;
			gsl_matrix_set(A,i,j,rnum);
			gsl_matrix_set(A,j,i,rnum); } 
	}
	gsl_matrix_memcpy(B,A);

	
	printf("Random symmetric matrix:\n"); printMatrix(A);

	gsl_matrix* evecs = gsl_matrix_calloc(dim,dim);
	gsl_vector* evals = gsl_vector_calloc(dim);
	int sweepCount = jacobi(A, evecs, evals);

	printf("Matrix after Jacobi diagonalization:\n");
	printMatrix(A); 
	printf("Number of sweeps: %i \n", sweepCount);
	printf("Corresponding eigenvalues:\n");
	gsl_vector_fprintf(stdout,evals,FORMAT);

	printf("Checking that V^T A V = D, where D is the matrix containing the eigenvalues on the diagonal:\n");
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,B,evecs,0,A);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,evecs,A,0,B);
	printMatrix(B);


	gsl_matrix_free(A); gsl_matrix_free(B); gsl_matrix_free(evecs); gsl_vector_free(evals);
	return 0;
}
