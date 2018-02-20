#include <stdio.h>
#include <tgmath.h>
#include <stdlib.h>
#include <math.h>
#include "function-precision.h"


int main(){
	int x;
	x=equal(1.0,3.0,0.0,1.0000000000001);
	printf("Return value= %i\n",x);
	return 0;
}

