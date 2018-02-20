#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char** argv)
{
	for(int i=1; i<argc; i++){
		double x = atof(argv[i]);
		printf("%lg \t %lg\n",x,sin(x));
	}
	return 0;
}