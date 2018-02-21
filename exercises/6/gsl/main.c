#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_airy.h>

int main () {
	//Plotting airy functions
	FILE* streamAi = fopen("outAi.txt","w");
	FILE* streamBi = fopen("outBi.txt","w");
	for(double i=-10; i<2; i=i+0.002){
		fprintf(streamAi, "%g \t %g\n",i, gsl_sf_airy_Ai(i,GSL_PREC_DOUBLE));
		fprintf(streamBi, "%g \t %g\n",i, gsl_sf_airy_Bi(i,GSL_PREC_DOUBLE));
	}
	fclose(streamAi);
	fclose(streamBi);



	/*Solving system of linear equations*/

	//Define matrix
	gsl_matrix* m=gsl_matrix_alloc(3,3);
	gsl_matrix_set(m,0,0, 6.13);
	gsl_matrix_set(m,0,1, -2.90);
	gsl_matrix_set(m,0,2, 5.86);
	gsl_matrix_set(m,1,0, 8.08);
	gsl_matrix_set(m,1,1, -6.31);
	gsl_matrix_set(m,1,2, -3.89);
	gsl_matrix_set(m,2,0, -4.36);
	gsl_matrix_set(m,2,1, 1.00);
	gsl_matrix_set(m,2,2, 0.19);

	printf("Elements of matrix m:\n");
	gsl_matrix_fprintf(stdout,m,"%g");
	printf("\n");

	//Define vector
	gsl_vector* v=gsl_vector_alloc(3);
	gsl_vector_set(v,0,6.23);
	gsl_vector_set(v,1,5.37);
	gsl_vector_set(v,2,2.29);

	printf("Elements of vector v:\n");
	gsl_vector_fprintf(stdout, v, "%g");
	printf("\n");

	gsl_vector* x = gsl_vector_alloc(3);
	gsl_linalg_HH_solve(m, v, x); //Solves m*x=v and stores result in x
	printf("Vector x s.t. m*x=v:\n");
	gsl_vector_fprintf(stdout,x,"%g");
	printf("\n");



	printf("m*x evaluates to:\n");
	gsl_matrix* m1=gsl_matrix_alloc(3,3);
	gsl_matrix_set(m1,0,0, 6.13);
	gsl_matrix_set(m1,0,1, -2.90);
	gsl_matrix_set(m1,0,2, 5.86);
	gsl_matrix_set(m1,1,0, 8.08);
	gsl_matrix_set(m1,1,1, -6.31);
	gsl_matrix_set(m1,1,2, -3.89);
	gsl_matrix_set(m1,2,0, -4.36);
	gsl_matrix_set(m1,2,1, 1.00);
	gsl_matrix_set(m1,2,2, 0.19);

	gsl_vector* v1 = gsl_vector_calloc(3);
	gsl_blas_dgemv(CblasNoTrans,1.0, m1, x, 0,v1);
	gsl_vector_fprintf(stdout,v1,"%g");
	printf("Succes!\n");

	gsl_matrix_free(m);
	gsl_matrix_free(m1);
	gsl_vector_free(x);
	gsl_vector_free(v);
	gsl_vector_free(v1);
	return 0;
}