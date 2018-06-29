#include <stdio.h>
#include <complex.h>
#include <tgmath.h>

int main() {
printf("First part of exercises:\n");
printf("gamma=%g\n",tgamma(5)); 
printf("bessel=%g\n", j1(0.5));
complex double x1=csqrt(-2);
printf("x1 = sqrt(-2) = %g + %gi\n", creal(x1), cimag(x1));
complex double x2=cexp(I);
printf("x2 = e^i = %g + %gi\n", creal(x2), cimag(x2));
complex double x3=cexp(I*M_PI);
printf("x3 = e^(i*pi) = %g + %gi\n", creal(x3), cimag(x3));
complex double x4=pow(I,M_E);
printf("x4 = i^e = %g + %gi\n", creal(x4), cimag(x4));

printf("Second part of exercises:\n");
float yfloat=0.111111111111111111111111111111;
printf("yfloat=%.25g\n",yfloat);
double ydouble=0.111111111111111111111111111111;
printf("ydouble=%.25lg\n", ydouble);
long double ylongdouble=0.111111111111111111111111111111L;
printf("ylongdouble=%.25Lg\n", ylongdouble);

return 0;
}

