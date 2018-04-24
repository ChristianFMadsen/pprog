#ifndef adapInt
#define adapInt
double int1(double x, int* fCalls);
double int2(double x, int* fCalls);
double int3(double x, int* fCalls);
double int4(double x, int* fCalls);
double int5(double x, int* fCalls);
double int6(double x, int* fCalls);
double intRoutine(double f(double x, int* fCalls),double a, double b, double acc, double eps, double x2, double x3, double* error, int numberOfRecursions, int* fCalls);
double my_integrator(double f(double x, int* fCalls),double a,double b, double acc,double eps, double* error, int* fCalls);
double aToInf(double f(double x, int* fCalls), double a, double acc, double eps, double* error, int* fCalls);
double clenCurtis(double f(double x, int* fCalls), double acc, double eps, double* error, int* fCalls);
#endif