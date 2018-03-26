#include <stdlib.h>
#include <assert.h>

typedef struct {int n; double *x, *y, *b, *c;} qspline;

//Function that builds quadratic spline
qspline* qspline_alloc(int n, double* x, double* y)
{
qspline* v=malloc(sizeof(qspline));


v->b = malloc((n-1)*sizeof(double));
v->c = malloc((n-1)*sizeof(double));
v->x = malloc(n*sizeof(double));
v->y = malloc(n*sizeof(double));
v->n = n;

for(int i=0; i<n; i++)
{
	v->x[i]=x[i];
	v->y[i]=y[i];
}

double p[n-1], dx[n-1], dy[n-1];

//Calculate p_i's:
for (int i = 0; i < n-1; i++)
{
	dx[i]=x[i+1]-x[i];
	dy[i]=y[i+1]-y[i];
	p[i]=dy[i]/dx[i];
}

//Assign c_1=0:
v->c[0]=0;

//Calculate c_i's recursively upwards:
for (int i = 0; i < n-2; i++)
{
	v->c[i+1] = (p[i+1]-p[i]-(v->c[i])*dx[i])/dx[i+1];
}

//Assign c_(n-1)=1/2*c_(n-1) using c_(n-1) from the upwards recursion:
v->c[n-2]=(v->c[n-2])/2.0;

//Calculate c_i's recursively downwards:
for (int i = n-3; i >= 0; i--)
{
	v->c[i]=(p[i+1]-p[i]-(v->c[i+1])*dx[i+1])/dx[i];
}

//Calculate b_i's:
for (int i = 0; i < n-1; i++)
{
	v->b[i]=p[i]-(v->c[i])*dx[i];
}


return v;
}

//Function that evaluates quadratic spline at a point z
double qspline_eval(qspline* v, double z)
{
//Assert point z is between x_inital and x_final:
assert(z>=v->x[0] && z<=v->x[v->n-1]);

//Binary search:
int i=0;
int j=v->n-1;
while(j-i>1)
{
	int m=(i+j)/2;
	if(z > v->x[m]) i=m;
	else j=m;
} 

double zxdiff=z-(v->x[i]);
return v->y[i]+zxdiff*(v->b[i]+zxdiff*v->c[i]);
}

//Function that free's the allocated memory
void qspline_free(qspline* v)
{
	free(v->x);
	free(v->y);
	free(v->b);
	free(v->c);
	free(v);
}

//Function to evaluate derivative at z
double qspline_derivative(qspline* v, double z)
{
assert(v->n>1);
assert(z>=v->x[0] && z<=v->x[v->n-1]); //Check that z is in the range of x-values

/*binary search*/
int i=0, j=v->n-1, k;

while(j-i>1)
{
k=(j+i)/2;
if(z>v->x[k]) i=k;
else j=k;
}

double deriv=v->b[i]+2.0*v->c[i]*(z-v->x[i]);
return deriv;
}

double qspline_integral(qspline* v, double z)
{
assert(v->n>1);
assert(z>=v->x[0] && z<=v->x[v->n-1]);

/*binary search*/
int i=0, j=v->n-1, k;

while(j-i>1)
{
k=(j+i)/2;
if(z>v->x[k]) i=k;
else j=k;
}

//Calculate integral of quadratic spline from x[0] to z
double ires;
double integral_result=0;
for (int m = 0; m <= i; m++)
{
	if(m<i)
	{
	ires=v->c[m]*v->x[m]*v->x[m]*v->x[m+1] + v->b[m]*v->x[m]*v->x[m]
	- v->c[m]*v->x[m]*v->x[m+1]*v->x[m+1] - v->b[m]*v->x[m]*v->x[m+1]
	- v->c[m]*v->x[m]*v->x[m]*v->x[m]/3.0 - v->b[m]*v->x[m]*v->x[m]/2.0
	- v->y[m]*v->x[m] + v->c[m]*v->x[m+1]*v->x[m+1]*v->x[m+1]/3.0
	+ v->b[m]*v->x[m+1]*v->x[m+1]/2.0 + v->y[m]*v->x[m+1];

	integral_result=integral_result+ires;
	}

	else if(m==i)
	{
	ires=v->c[m]*v->x[m]*v->x[m]*z + v->b[m]*v->x[m]*v->x[m]
	- v->c[m]*v->x[m]*z*z - v->b[m]*v->x[m]*z
	- v->c[m]*v->x[m]*v->x[m]*v->x[m]/3.0 - v->b[m]*v->x[m]*v->x[m]/2.0
	- v->y[m]*v->x[m] + v->c[m]*z*z*z/3.0
	+ v->b[m]*z*z/2.0 + v->y[m]*z;

	integral_result=integral_result+ires;
	}

}


return integral_result;
}