#include <stdio.h>
#include <gsl/gsl_matrix.h>
#define FORMAT "%8.3f"


void printMatrix(gsl_matrix* A){
	for(int i=0;i<A->size1;i++){
		for(int j=0;j<A->size2;j++) printf(FORMAT,gsl_matrix_get(A,i,j));
		printf("\n"); }
}