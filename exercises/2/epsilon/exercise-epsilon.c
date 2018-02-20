#include <stdio.h>
#include <limits.h>
#include <float.h>

int main()
{
	
	//Exercise 1:
	printf("Exercise 1:\n");
	// INT_MIN stuff:
	int x1=INT_MAX;
	printf("Max integer according to limits.h: INT_MAX=%i\n",x1);

	int i=1;
	while(i+1>i) {
		i++;
	}
	printf("My INT_MAX_WHILE=%i\n",i);

	int j;
	for (j=1; j+1 > j; ++j){

	}
	printf("My INT_MAX_FOR=%i\n", j);

	int k=1;
	do{
		k++;
	} while(k+1>k);
	printf("My INT_MAX_DO_WHILE=%i\n", k);


	//INT_MIN stuff:
	int x2=INT_MIN;
	printf("Min integer according to limits.h: INT_MIN=%i\n",x2);

	int l=1;
	while(l-1<l) {
		l--;
	}
	printf("My INT_MIN_WHILE=%i\n",l);

	int h;
	for (h=1; h-1 < h; h--){

	}
	printf("My INT_MIN_FOR=%i\n", h);

	int g=1;
	do{
		g--;
	} while(g-1<g);
	printf("My INT_MIN_DO_WHILE=%i\n", g);


	//Machine epsilon stuff:
	double e1=DBL_EPSILON;
	printf("Machine epsilon as defined in float.h (double): epsilon=%g\n", e1);

	double e2=1.0;
	while(1.0+e2!=1.0){
		e2 = e2/2;
	}
	printf("my_epsilon_dbl=%g\n",2*e2);

	float e3=FLT_EPSILON;
	printf("Machine epsilon as defined in float.h (float): epsilon=%g\n", e3);

	float e4=1.0f;
	while(1.0f+e4!=1.0f){
		e4 = e4/2;
	}
	printf("my_epsilon_flt=%g\n",2*e4);

	long double e5 = LDBL_EPSILON;
	printf("Machine epsilon as defined in float.h (lng dbl): epsilon=%Lg\n", e5);

	long double e6=1.0L;
	while(1.0L+e6!=1.0L){
		e6=e6/2;
	}
	printf("my_epsilon_lngdbl=%Lg\n",2*e6);

	// Exercise 2:
	printf("Exercise 2:\n");
	
	int max=INT_MAX/2;
	int c_up=1;
	float sum_up_float=0.0f;
	while(c_up <= max){
		sum_up_float=sum_up_float +1.0f/c_up;
		c_up++;
	}
	printf("sum_up_float=%g\n", sum_up_float);

	int c_down=max;
	float sum_down_float = 0.0f;
	while(c_down >= 1){
		sum_down_float=sum_down_float + 1.0f/c_down;
		c_down--;
	}
	printf("sum_down_float=%g\n", sum_down_float);
	printf("difference_float=%g\n", sum_down_float-sum_up_float);

	printf("Explain the difference:\n");
	printf("The floating point errors are bigger for bigger numbers so sum_up_float has a bigger error when calculating the sum because it starts with the biggest numbers which means bigger errors start to accumulate earlier.\n");

	printf("Does this sum converge as a function of max?\n");
	printf("No. Making max bigger does not reduce difference_float but making max significantly smaller does. This suggests that the sum diverges as max -> INT_MAX.\n");

	int c_upd=1;
	double sum_up_double=0.0;
	while(c_upd <= max){
		sum_up_double=sum_up_double +1.0/c_upd;
		c_upd++;
	}
	printf("sum_up_double=%.25g\n", sum_up_double);

	int c_downd=max;
	double sum_down_double = 0.0;
	while(c_downd >= 1){
		sum_down_double=sum_down_double + 1.0/c_downd;
		c_downd--;
	}
	printf("sum_down_double=%.25g\n", sum_down_double);
	printf("difference_double=%g\n", sum_up_double-sum_down_double);

	printf("Explain the result:\n");
	printf("The double data type is more precise i.e. it contains more digits which means the sums are almost identical. There is still an error but it is very small compared to before. The error can be explained in the same way as before.\n");



	return 0;
}
