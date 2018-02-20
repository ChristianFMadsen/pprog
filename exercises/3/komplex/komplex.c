#include <stdio.h>
#include "komplex.h"
#include <math.h>
#include <complex.h>

void komplex_print(char* s, komplex z){
	printf("%s (%g,%g)\n", s, z.re, z.im);
}

void komplex_set(komplex* z, double x, double y){
	/* printf("(%g,%g) has been redefined to:\n", z.re, z.im); */
	(*z).re = x; /*Redefine real part*/
	(*z).im = y; /*Redefine imaginary part*/
	/* printf("(%g,%g)\n", z.re, z.im); */
}

komplex komplex_new (double x, double y){
	komplex z={x,y}; /*Change value of real part to x and value of imaginary part to y*/
	return z;
}

komplex komplex_add(komplex x, komplex y){
	komplex result = {x.re+y.re,x.im+y.im};
	return result;
}

komplex komplex_sub(komplex x, komplex y){
	komplex result = {x.re-y.re,x.im-y.im};
	return result;
}

int komplex_equal(komplex x, komplex y){
	if (x.re==y.re && x.im==y.im){
		return 1;
	}

	else{
		return 0;
	}
}

komplex komplex_mul(komplex x, komplex y){
	komplex result = {x.re*y.re-x.im*y.im, x.im*y.re+x.re*y.im};
	return result;
}

komplex komplex_div(komplex x, komplex y){
	komplex result = {(x.re*y.re+x.im*y.im)/(pow(y.re,2)+pow(y.im,2)),(x.im*y.re-x.re*y.im)/(pow(y.re,2)+pow(y.im,2))};
	return result;
}

komplex komplex_conjugate(komplex z){
	komplex result = {z.re, -z.im};
	return result;
}

komplex komplex_abs(komplex z){
	komplex result = {pow( pow(z.re,2) + pow(z.im,2),0.5),0};
	return result;
}

komplex komplex_exp(komplex z){
	komplex result = {exp(z.re)*cos(z.im),exp(z.re)*sin(z.im)};
	return result;
}

komplex komplex_cos(komplex z){
	komplex result = {creal(1.0/2*(cexp(I*(z.re+I*z.im))+cexp(-I*(z.re+I*z.im)))), cimag(1.0/2*(cexp(I*(z.re+I*z.im))+cexp(-I*(z.re+I*z.im))))};
	return result;
}

komplex komplex_sin(komplex z){
	komplex result = {creal(1.0/(2*I)*(cexp(I*(z.re+I*z.im))-cexp(-I*(z.re+I*z.im)))), cimag(1.0/(2*I)*(cexp(I*(z.re+I*z.im))-cexp(-I*(z.re+I*z.im))))};
	return result;
}

komplex komplex_sqrt(komplex z){
	double r;
	double theta;
	r=sqrt(pow(z.re,2)+pow(z.im,2));
	theta=atan2(z.im, z.re);
	komplex result = {creal(sqrt(r)*cexp(I*theta/2)), cimag(sqrt(r)*cexp(I*theta/2))};
	return result;
}