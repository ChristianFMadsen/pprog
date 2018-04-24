#ifndef ode
#define ode
void rkstep12(double x, double h, gsl_vector* yx, void f(double x, gsl_vector* y, gsl_vector* dydx), gsl_vector* yxh, gsl_vector* err);
int ode_driver(gsl_vector* xlist, gsl_matrix* ylist, double b, double h, double acc, double eps, int maxSteps, int dim, void f(double x, gsl_vector* y,gsl_vector* dydx));
void sinODE(double x, gsl_vector* y, gsl_vector* dydx);
double integrator(gsl_vector* xlist, gsl_matrix* ylist, double b, double h, double acc, double eps, int maxSteps, int dim, double f(double x));
double intFun(double x);
#endif