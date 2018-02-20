#include "nvector.h"
#include <stdio.h>
#include <stdlib.h>
#define RND (double) rand()/RAND_MAX

int main(){
	int n=7;

	printf("Testing nvector_alloc:\n");
	nvector *v = nvector_alloc(n);
	if (v == NULL){
		printf("Test failed\n");
	}
	else{
		printf("Test passed\n");
	}

	printf("Testing nvector_set and nvector_get:\n");
	double value = RND; /*Pseudo random number between 0 and 1*/
	int i = n/2; /* If n=7 this returns i=3? */
	nvector_set(v, i, value);
	double vi = nvector_get(v,i);
	if(vi==value){
		printf("Test passed\n");
	}
	else{
		printf("Test failed\n");
	}


	printf("Testing nvector_dot_product:\n");
	nvector *a = nvector_alloc(n);
	nvector *b = nvector_alloc(n);
	nvector *d = nvector_alloc(n);
	nvector *e = nvector_alloc(n);
	nvector *e_2 = nvector_alloc(n); //Used in nvector_scale. 2*e
	nvector *ab_add = nvector_alloc(n); //Used in nvector_add
	nvector *de_sub = nvector_alloc(n); //Used in nvector_sub
	nvector *c = nvector_alloc(n+1);
	/*Fills vectors a and b with random numbers*/
	for (int i=0; i<=n-1; i++){
		double x=RND;
		double y=RND;
		nvector_set(a, i, x);
		nvector_set(b, i, y);
		nvector_set(d, i, x);
		nvector_set(e, i, y);
		nvector_set(ab_add, i, x+y);
		nvector_set(de_sub, i, x-y);
		nvector_set(e_2, i, 2*y);
	}

	/*Fills vector c with random numbers*/
	for (int i=0; i<=n; i++){
		double z=RND;
		nvector_set(c, i, z);
	}

	printf("Dot product between vectors of same dimension:\n");
	nvector_dot_product(a,b);
	printf("Dot product between vectors of different dimensions:\n");
	nvector_dot_product(a,c);

	printf("Testing nvector_print:\n");
	nvector_print("Contents of vector a", a);
	nvector_print("Contents of vector b",b);
	nvector_print("Contents of vector c",c);
	printf("Succes!\n");

	printf("Testing nvector_set_zero:\n");
	nvector_set_zero(c);
	nvector_print("nvector_set_zero(c) returns:",c);
	printf("Succes!\n");

	printf("Testing nvector_equal:\n");
	printf("nvector_equal will return 1 if the vectors are equal and 0 otherwise.\n");
	printf("nvector_equal(a,a) returns:\n");
	printf("%i\n", nvector_equal(a,a));
	printf("nvector_equal(a,b) returns:\n");
	printf("%i\n", nvector_equal(a,b));
	printf("nvector_equal(a,c) returns:\n");
	printf("%i\n", nvector_equal(a,c));
	printf("Succes!\n");

	printf("Testing nvector_add:\n");
	nvector_print("Contents of vector a", a);
	nvector_print("Contents of vector b",b);
	nvector_print("a+b is:", ab_add);
	nvector_add(a,b);
	nvector_print("nvector_add(a,b) returns:", a);
	printf("nvector_add(a,c) returns:\n");
	nvector_add(a,c);
	printf("Succes!\n");

	printf("Testing nvector_sub:\n");
	nvector_print("Contents of vector d", d);
	nvector_print("Contents of vector e",e);
	nvector_print("d-e is:", de_sub);
	nvector_sub(d,e);
	nvector_print("nvector_sub(d,e) returns:", d);
	printf("nvector_sub(d,c) returns:\n");
	nvector_sub(d,c);
	printf("Succes!\n");

	printf("Testing nvector_scale:\n");
	nvector_print("Contents of vector e",e);
	nvector_print("2*e is:", e_2);
	nvector_scale(e,2);
	nvector_print("nvector_scale(e,2) returns:",e);
	printf("Succes!\n");



return 0;
}