#include "komplex.h"
#include <stdio.h>
#include <math.h>
#include <complex.h>

int main(){

	komplex a={1,2};
	komplex b={3,4};
	
	komplex_print("a =",a);
	komplex_print("b =",b);
	
	/*
	komplex c={10,-3};
	komplex_print("c =",c);
	komplex_set(c,1,1);
	*/


	printf("Testing komplex_add:\n");
	komplex r1 = komplex_add(a,b);
	komplex R1 = {4,6};
	komplex_print("a+b should equal:",R1);
	komplex_print("komplex_add returns:",r1);
	printf("Succes!\n");

	printf("Testing komplex_sub:\n");
	komplex r2 = komplex_sub(a,b);
	komplex R2 = {-2,-2};
	komplex_print("a-b should equal:",R2);
	komplex_print("komplex_sub returns:",r2);
	printf("Succes!\n");

	printf("Testing komplex_equal:\n");
	int x2 = komplex_equal(a,a);
	int x3 = komplex_equal(a,b);
	printf("Return value for komplex_equal(a,a): %i\n",x2);
	printf("Return value for komplex_equal(a,b): %i\n",x3);
	printf("Succes!\n");

	printf("Testing komplex_new:\n");
	komplex n = komplex_new(1.2,-3.9);
	komplex_print("Input komplex_new(1.2,-3.9). Output:", n);
	printf("Succes!\n");

	printf("Testing komplex_mul:\n");
	komplex r3 = komplex_mul(a,b);
	komplex R3 = {-5,10};
	komplex_print("a*b should equal:",R3);
	komplex_print("komplex_mul returns:", r3);
	printf("Succes!\n");

	printf("Testing komplex_div:\n");
	komplex r4 = komplex_div(a,b);
	komplex R4 = {11.0/25.0,2.0/25.0};
	komplex_print("a/b should equal:",R4);
	komplex_print("komplex_div returns:", r4);
	printf("Succes!\n");

	printf("Testing komplex_conjugate:\n");
	komplex r5 = komplex_conjugate(a);
	komplex R5 = {a.re,-a.im};
	komplex_print("The complex conjugate of a is:", R5);
	komplex_print("komplex_conjugate returns:", r5);
	printf("Succes!\n");

	printf("Testing komplex_abs:\n");
	komplex r6 = komplex_abs(a);
	komplex R6 = {pow( pow(a.re,2) + pow(a.im,2) ,0.5),0};
	komplex_print("The absolute value of the complex number a is:",R6);
	komplex_print("komplex_abs returns:",r6);
	printf("Succes!\n");

	printf("Testing komplex_exp:\n");
	komplex r7 = komplex_exp(a);
	komplex R7 = {exp(a.re)*cos(a.im),exp(a.re)*sin(a.im)};
	komplex_print("e^(1+2i) is equal to:",R7);
	komplex_print("komplex_exp returns:",r7);
	printf("Succes!\n");

	printf("Testing komplex_cos:\n");
	komplex r8 = komplex_cos(a);
	komplex R8 = {creal(ccos(a.re+a.im*I)),cimag(ccos(a.re+a.im*I))};
	komplex_print("cos(1+2i) is equal to:",R8);
	komplex_print("komplex_cos returns:",r8);
	printf("Succes!\n");

	printf("Testing komplex_sin:\n");
	komplex r9 = komplex_sin(a);
	komplex R9 = {creal(csin(a.re+a.im*I)), cimag(csin(a.re+a.im*I))};
	komplex_print("sin(1+2i) is equal to:", R9);
	komplex_print("komplex_sin returns:", r9);

	printf("Testing komplex_sqrt:\n");
	komplex r10 = komplex_sqrt(a);
	komplex R10 = {creal(csqrt(a.re+a.im*I)), cimag(csqrt(a.re+a.im*I))};
	komplex_print("sqrt(1+2i) is equal to:", R10);
	komplex_print("komplex_sqrt returns:", r10);
	printf("Succes!\n");


}