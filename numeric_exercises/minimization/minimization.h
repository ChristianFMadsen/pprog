#ifndef minimization
#define minimization
int newton(double f(gsl_vector* x, gsl_vector* df, gsl_matrix* H), gsl_vector* xstart, double eps);
int qNewton(double f(gsl_vector* x), gsl_vector* xstart, double eps, double dx);
void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R);
void qr_gs_solve(const gsl_matrix* Q, const gsl_matrix* R, const gsl_vector* b, gsl_vector* x);
void qr_gs_inverse(const gsl_matrix* Q, const gsl_matrix* R, gsl_matrix* B);
double f_ros(gsl_vector* x1, gsl_vector* df, gsl_matrix* H);
double f_him(gsl_vector* x1, gsl_vector* df, gsl_matrix* H);
double fx_him(gsl_vector* x1);
double fx_ros(gsl_vector* x1);
void gradient(double f(gsl_vector* x), gsl_vector* x, gsl_vector* df, double dx);
double decay_fun(gsl_vector* x);
#endif