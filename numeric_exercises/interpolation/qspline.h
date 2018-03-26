#ifndef QSPLINE_H
#define QSPLINE_H
typedef struct {int n; double *x, *y, *b, *c;} qspline;
qspline* qspline_alloc(int n, double* x, double* y);
double qspline_eval(qspline* v, double z);
double qspline_derivative(qspline* v, double z); 
double qspline_integral(qspline* v, double z);  
void qspline_free(qspline* v);
#endif
