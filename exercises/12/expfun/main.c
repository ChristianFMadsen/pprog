#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

double myexp(double x);

int main()
{
	double x_start=-5;
	double x_end=5;
	double dx=0.1;

	printf("#x \t\t myexp \t\t math.h exp\n");
	for (double j=x_start; j<x_end+dx; j=j+dx)
	{
		printf("%.6e \t %.6e \t %.6e\n",j, myexp(j), exp(j));
	}
	return 0;
}