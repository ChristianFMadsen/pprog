#ifndef QR_LS
#define QR_LS
void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R);
void qr_gs_solve(const gsl_matrix* Q, const gsl_matrix* R, const gsl_vector* b, gsl_vector* x);
void qr_gs_inverse(const gsl_matrix* Q, const gsl_matrix* R, gsl_matrix* B);
double funs(int i, double x);
void leastsquares(double funs(int i, double x), int colDim, gsl_vector* x,  gsl_vector* y,  gsl_vector* errors,  gsl_vector* c,  gsl_matrix* S);
#endif
