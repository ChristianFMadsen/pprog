#include <stdio.h>
#include <tgmath.h>
#include <stdlib.h>
#include <math.h>
#include "function-precision.h"

int equal(double a, double b, double tau, double epsilon){
	if (fabs(a-b)<tau) {
		//printf("absolute\n");
		return 1;
	}
	else if( ( fabs(a-b) ) / ( fabs(a) + fabs(b) ) < epsilon/2){
		//printf("relative\n");
		return 1;
	}
	else{
		//printf("else\n");
	return 0;
	}	
}