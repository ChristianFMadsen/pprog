#ifndef mc
#define mc
double int1(double* x);
double int2(double* x);
double int3(double* x);
void rndx(int dim, double *a, double *b, double *x);
void plainMC(int dim, double *a, double *b, double f(double* x), int N, double* result, double* error);
#endif