#include<math.h>
#include<assert.h>
#include"montecarlo.h"


double intRoutine(double f(double x, int* fCalls),double a, double b, double acc, double eps, double x2, double x3, double* error, int numberOfRecursions, int* fCalls){
	assert(numberOfRecursions<1e6); //Stop if the integral still hasn't converged after a million recursions.
	double x1=f(a+(b-a)/6,fCalls), x4=f(a+5*(b-a)/6,fCalls); //x1 and x4 calculated according to eq. 48
	double Q=(2*x1+x2+x3+2*x4)/6*(b-a); //Q defined in eq. 44 with the weights defined in eq. 49
	double q=(x1+x4+x2+x3)/4*(b-a); //q defined in eq. 45, weights defined in eq. 50
	double tolerance=acc+eps*fabs(Q); //Calculate the tolerance according to eq. 47 
	*error=fabs(Q-q); //calculate error estimate according to eq. 46  
	if(*error < tolerance){
		return Q; } //If error<tolerance the integral has been calculated satisfactorily -> return the integral value. 
	else {
		double Q1=intRoutine(f,a,(a+b)/2,acc/sqrt(2.0),eps,x1,x2,error,numberOfRecursions+1,fCalls); //divide into smaller interval to see if it converges, and adjust accuracy by 1/sqrt(2) as found under eq. 47
		double Q2=intRoutine(f,(a+b)/2,b,acc/sqrt(2.0),eps,x3,x4,error,numberOfRecursions+1,fCalls); //same idea as above just from the other side (right of the middle)
		return Q1+Q2;
	}
}

double my_integrator(double f(double x, int* fCalls),double a,double b, double acc,double eps, double* error, int* fCalls){
	double x2=f(a+2*(b-a)/6,fCalls), x3=f(a+4*(b-a)/6,fCalls); //x2 and x3 calculated according to eq. 48
	int numberOfRecursions=0; //reset number of recursions
	return intRoutine(f,a,b,acc,eps,x2,x3,error,numberOfRecursions,fCalls);
}



double aToInf(double f(double x, int* fCalls), double a, double acc, double eps, double* error, int* fCalls){
	double g(double t, int* fCalls){
		return f((a + (1-t)/t), fCalls)*pow(t,-2);
	}
	return my_integrator(g,0,1,acc,eps,error,fCalls);
}


double clenCurtis(double f(double x, int* fCalls), double a, double b, double acc, double eps, double* error, int* fCalls){
	double g(double t, int* fCalls){
		return f((a+b)/2+(a-b)/2*cos(t), fCalls)*sin(t)*(b-a)/2;
	}
	return my_integrator(g,0,M_PI,acc,eps,error,fCalls);
}