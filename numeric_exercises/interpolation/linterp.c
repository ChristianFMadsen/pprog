#include <assert.h>

double linterp(int n, double* x, double* y, double z)
{
assert(n>1);
assert(z>=x[0] && z<=x[n-1]); //Check that z is in the range of x-values

/*binary search*/
int i=0, j=n-1, k;

while(j-i>1)
{
k=(j+i)/2;
if(z>x[k]) i=k;
else j=k;
}

/*calculate linear spline*/
return y[i]+(y[i+1]-y[i])/(x[i+1]-x[i])*(z-x[i]);
} 


double linterp_integ(int n, double *x, double *y, double z)
{
assert(n>1);
assert(z>=x[0] && z<=x[n-1]); //Check that z is in the range of x-values


/*binary search*/
int i=0, j=n-1, k;

while(j-i>1)
{
k=(j+i)/2;
if(z>x[k]) i=k;
else j=k;
}

//Calculate integral of linear spline from x[0] to z
double pi, ires;
double integral_result=0;
for (int c = 0; c <= i; c++)
{
	if(c<i)
	{
	pi=(y[c+1]-y[c])/(x[c+1]-x[c]);
	ires=(x[c+1]-x[c])*(2*y[c]+pi*x[c]-2*pi*x[c]+pi*x[c+1])/2.0;
	integral_result=integral_result+ires;
	}

	else if(c==i)
	{
		pi=(y[c+1]-y[c])/(x[c+1]-x[c]);
		ires=(z-x[c])*(2*y[c]+pi*x[c]-2*pi*x[c]+pi*z)/2.0;
		integral_result=integral_result+ires;
	}

}


return integral_result;
}