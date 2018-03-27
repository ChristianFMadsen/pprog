#ifndef LIN_EQ
#define LIN_EQ
void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R);
void qr_gs_solve(const gsl_matrix* Q, const gsl_matrix* R, const gsl_vector* b, gsl_vector* x);
void qr_gs_inverse(const gsl_matrix* Q, const gsl_matrix* R, gsl_matrix* B);
#endif
