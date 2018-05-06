#include <gsl/gsl_matrix.h>
#ifndef have_ann
#define have_ann
typedef struct {int n; double (*f)(double x); gsl_vector* data;} ann;
ann* ann_alloc(int n, double (*f)(double x)); 
void ann_free(ann* network);
double ann_feed_forward(ann* network, double x);
void ann_train(ann* network, gsl_vector* vx, gsl_vector* vf);
int qNewton(double f(gsl_vector* x), gsl_vector* xstart, double eps, double dx);
void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R);
void qr_gs_solve(const gsl_matrix* Q, const gsl_matrix* R, const gsl_vector* b, gsl_vector* x);
void qr_gs_inverse(const gsl_matrix* Q, const gsl_matrix* R, gsl_matrix* B);
void gradient(double f(gsl_vector* x), gsl_vector* x, gsl_vector* df, double dx);
#endif