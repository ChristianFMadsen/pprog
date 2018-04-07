#ifndef md
#define md
int jacobi(gsl_matrix* A, gsl_matrix* evecs, gsl_vector* evals);
void printMatrix(gsl_matrix* A);
int jacobi_eig_descending(gsl_matrix* A, gsl_vector* evalsVec, int k);
int jacobi_eig_ascending(gsl_matrix* A, gsl_vector* evalsVec, int k);
#endif