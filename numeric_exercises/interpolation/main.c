#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "lspline.h"
#include "qspline.h"

int main() 
{
FILE* dataStream=fopen("data.txt","w");

int n1=600;
double x[n1];
double y[n1];
fprintf(dataStream,"#x \t y \t y_int \t y_deriv\n");
for (int i = 0; i < n1; i++)
{
	x[i]=i/50.0;
	y[i]=sin(5*x[i]);
	fprintf(dataStream,"%g \t %g \t %g \t %g\n", x[i], y[i], 
		2/5.0*sin(5*x[i]/2.0)*sin(5*x[i]/2.0), 5.0*cos(5.0*x[i]));
}

fprintf(dataStream,"\n\n");
qspline* q=qspline_alloc(n1,x,y);
fprintf(dataStream,"#x \t y (lspline) \t y (qspline) \t gsl_int \t lspline_int \t qspline_int \t qspline_deriv\n");
for (double z=x[0]; z<=x[n1-1]; z=z+0.01)
{
	double qspl=qspline_eval(q,z);
	fprintf(dataStream,"%g \t %g \t %g \t %g \t %g \t %g \t %g\n",z, linterp(n1,x,y,z), qspl,
	sinintegral(0.0,z), linterp_integ(n1, x, y, z), qspline_integral(q, z), qspline_derivative(q, z));
}

printf("Comparison between manually computed quadratic spline coefficients b_i and c_i, and qspline.c:\n");
int n2=5;
double x1[n2], y1[n2], y2[n2], y3[n2];
x1[0]=1; x1[1]=2; x1[2]=3; x1[3]=4; x1[4]=5; 
y1[0]=y1[1]=y1[2]=y1[3]=y1[4]=1;
y2[0]=1; y2[1]=2; y2[2]=3; y2[3]=4; y2[4]=5; 
y3[0]=1; y3[1]=4; y3[2]=9; y3[3]=16; y3[4]=25;

qspline* q1=qspline_alloc(n2, x1, y1);
qspline* q2=qspline_alloc(n2, x1, y2);
qspline* q3=qspline_alloc(n2, x1, y3);

printf("Expected coefficients b_i and c_i:\n");
printf("b_11=b_12=b_13=b_14=c_11=c_12=c_13=c_14=0\n");
printf("b_21=b_22=b_23=b_24=1 & c_21=c_22=c_23=c_24=0\n");
printf("b_31=2, b_32=4, b_33=6, b_34=8 & c_31=c_32=c_33=c_34=1\n");

printf("Computed coefficients:\n");
for (int i = 0; i < n2-1; i++)
{
	printf("i=%i: b_1i=%g, b_2i=%g, b_3i=%g, c_1i=%g, c_2i=%g, c_3i=%g\n", 
		i+1, q1->b[i], q2->b[i], q3->b[i], q1->c[i], q2->c[i], q3->c[i]);
}

qspline_free(q);
qspline_free(q1);
qspline_free(q2);
qspline_free(q3);
fclose(dataStream);
return 0;
}