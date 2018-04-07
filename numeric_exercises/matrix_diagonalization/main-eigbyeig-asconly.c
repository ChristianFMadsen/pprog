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
	int numberOfEigVals=atoi(argv[2]);
	printf("Matrix dimension=%i\n", dim);

	gsl_matrix* A = gsl_matrix_calloc(dim,dim);
	gsl_vector* evals_asc = gsl_vector_calloc(numberOfEigVals);

	for (int i = 0; i < dim; i++){
		for (int j = i; j < dim; j++){ 
			double rnum = RND;
			gsl_matrix_set(A,i,j,rnum);
			gsl_matrix_set(A,j,i,rnum); } 
	}

	printf("Random symmetric matrix:\n"); printMatrix(A);
	printf("Matrix after having found the %i lowest eigenvalues:\n", numberOfEigVals);
	jacobi_eig_ascending(A,evals_asc,numberOfEigVals); printMatrix(A);
	printf("Found eigenvalues in ascending order (0's correspond to eigenvalues that haven't been calculated):\n");
	gsl_vector_fprintf(stdout, evals_asc, FORMAT);
	printf("Note that these are the correct eigenvalues by comparing to part A of the exercise.\n");

	gsl_matrix_free(A); gsl_vector_free(evals_asc);
	return 0;
}