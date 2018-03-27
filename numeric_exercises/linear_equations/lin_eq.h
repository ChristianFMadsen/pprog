#ifndef LIN_EQ
#define LIN_EQ
void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R);
void qr_gs_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
void qr_gs_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);
#endif
