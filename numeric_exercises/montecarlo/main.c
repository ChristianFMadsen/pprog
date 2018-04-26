#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "montecarlo.h"
#define RND ((double)rand()/RANDMAX)


int main(){
printf("Testing plain Monte Carlo integrator:\n");
int dim=1, N=1e6;
double a[dim], b[dim], result, error;
a[0]=-M_PI/2; b[0]=M_PI/2; 
printf("Testing on the integral tan(x) from %g to %g with N=%i. The result should be 0.\n", a[0],b[0],N);
printf("(actually the integral is not defined, however if you take the integral from -a to a and let a approach pi/2 then that limit is equal to 0.\n");
plainMC(dim,a,b,int1,N,&result,&error);
printf("integral_result=%g, estimated error=%g\n\n", result,error);
a[0]=0; b[0]=1;
printf("Now testing on the integral exp(x^2) from %g to %g with N=%i. The result should be approximately 1.46265...\n",a[0],b[0],N);
plainMC(dim,a,b,int2,N,&result,&error);
printf("integral_result=%g, estimated error=%g\n\n", result,error);
printf("Now testing on the ''difficult singular integral'' from the exercise page with N=%i. \n",N);
dim=3;
double a1[dim], b1[dim];
a1[0]=a1[1]=a1[2]=0;
b1[0]=b1[1]=b1[2]=M_PI;
plainMC(dim,a1,b1,int3,N,&result,&error);
printf("integral_result=%.10g, estimated error=%g\n", result,error);

dim=1;
FILE* datastream=fopen("errors.txt","w");
fprintf(datastream, "#N \t estimated error\n");
for (int i=1; i < 1000; i++){
	N=10*i;
	plainMC(dim,a,b,int2,N,&result,&error);
	fprintf(datastream, "%i \t %g\n",N,error);
}

fclose(datastream);
return 0;
}