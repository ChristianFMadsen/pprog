#include <gsl/gsl_matrix.h>
#ifndef have_ann
#define have_ann
typedef struct {int n; double (*g)(double); double (*dg)(double); gsl_vector* data;} ann;
ann* ann_alloc(int n, double (*g)(double), double (*dg)(double));
void ann_free(ann* network);
void ann_feed_forward(ann* network, double x, double* F, double* dF);
void ann_train(ann* network, gsl_vector* vx, double (*func)(double, double), double x0, double y0);
int qNewton(double f(gsl_vector* x), gsl_vector* xstart, double eps, double dx);
void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R);
void qr_gs_inverse(const gsl_matrix* Q, const gsl_matrix* R, gsl_matrix* B);
#endif