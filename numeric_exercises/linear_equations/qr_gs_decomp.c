#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "lin_eq.h"


void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R)
{
//A->Q (ie. A is destroyed) and R is build
int m = A->size2; //Number of columns in A
for (int i = 0; i < m; i++)
	{
		gsl_vector_view i_col = gsl_matrix_column(A,i); //i'th column
		double col_norm = gsl_blas_dnrm2(&i_col.vector); //euclidean norm of i'th column
		gsl_matrix_set(R,i,i,col_norm); //R_ii=sqrt(a_i^(TRANSPOSE)*a_i)
		gsl_vector_scale(&i_col.vector, 1.0/col_norm); //Normalization of column

		for (int j = i+1; j < m; j++)
			{
				gsl_vector_view j_col = gsl_matrix_column(A,j);
				double dot_result=0;
				gsl_blas_ddot(&i_col.vector, &j_col.vector, &dot_result);
				gsl_blas_daxpy(-dot_result, &i_col.vector, &j_col.vector); //j_col -> -dot_result*i_col + j_col (p. 3 just under GS ortho)
				gsl_matrix_set(R,i,j,dot_result);
			}	

	}	
}

void qr_gs_solve(const gsl_matrix* Q, const gsl_matrix* R, const gsl_vector* b, gsl_vector* x)
{
	int m = R->size1; //Number of rows in R
	gsl_blas_dgemv(CblasTrans, 1.0, Q, b, 0.0, x); //Q^(TRANSPOSED)*b stored in x
	/*Back substitution for loop */
	for (int i = m-1; i >= 0; i--)
	{
		double p=0;
		for (int k = i+1; k < m; k++)
		{
			p = p + gsl_matrix_get(R, i, k) * gsl_vector_get(x, k);
		}
		gsl_vector_set(x, i, (gsl_vector_get(x, i)-p)/gsl_matrix_get(R, i, i));
	}
}

void qr_gs_inverse(const gsl_matrix* Q, const gsl_matrix* R, gsl_matrix* B)
{
int m = Q->size1;
gsl_vector* b = gsl_vector_calloc(m);
gsl_vector* x = gsl_vector_calloc(m);
for (int i = 0; i < m; i++)
	{
		gsl_vector_set(b, i, 1.0);
		qr_gs_solve(Q, R, b, x);
		gsl_matrix_set_col(B, i, x);
		gsl_vector_set(b, i, 0.0);
	}
gsl_vector_free(b); gsl_vector_free(x);	
}
