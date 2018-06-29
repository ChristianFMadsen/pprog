#ifndef mc
#define mc
double int1(double* x);
double int2(double* x);
double int3(double* x);
void rndx(int dim, double *a, double *b, double *x);
void plainMC(int dim, double *a, double *b, double f(double* x), int N, double* result, double* error);
double intRoutine(double f(double x, int* fCalls),double a, double b, double acc, double eps, double x2, double x3, double* error, int numberOfRecursions, int* fCalls);
double my_integrator(double f(double x, int* fCalls),double a,double b, double acc,double eps, double* error, int* fCalls);
double aToInf(double f(double x, int* fCalls), double a, double acc, double eps, double* error, int* fCalls);
double clenCurtis(double f(double x, int* fCalls), double a, double b, double acc, double eps, double* error, int* fCalls);
double int2D(double f(double* x, int* fCalls), int* fCalls, double a, double b, double c(double), double d(double), double acc, double eps, double* error);
#endif