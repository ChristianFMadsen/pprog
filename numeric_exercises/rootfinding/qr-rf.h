#ifndef QR_RF
#define QR_RF
void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R);
void qr_gs_solve(const gsl_matrix* Q, const gsl_matrix* R, const gsl_vector* b, gsl_vector* x);
void qr_gs_inverse(const gsl_matrix* Q, const gsl_matrix* R, gsl_matrix* B);
int newton(int f(gsl_vector* x1, gsl_vector* fx), gsl_vector* x, double dx, double eps);
int newton_with_jacobian(int f(const gsl_vector* x, gsl_vector* fx, gsl_matrix* J), gsl_vector* x, double eps);
int newton_with_jacobian_refined(int f(const gsl_vector* x, gsl_vector* fx, gsl_matrix* J), gsl_vector* x, double eps);
#endif