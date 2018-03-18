#include <stdio.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

struct experimental_data {int n; double* t; double* y; double* e;};

double min_fun(const gsl_vector* x, void* params)
{
	double  A = gsl_vector_get(x,0);
	double  T = gsl_vector_get(x,1);
	double  B = gsl_vector_get(x,2);
	struct experimental_data* p = (struct experimental_data*) params;
	int     n = p->n;
	double* t = p->t;
	double* y = p->y;
	double* e = p->e;
	double sum=0;
	#define f(t) A*exp(-(t)/T) + B
	for(int i=0;i<n;i++)
	{
	 sum = sum + pow((f(t[i]) - y[i])/e[i],2);
	}
	return sum;
}


int main()
{
	//Data:
	double t[]= {0.47,1.41,2.36,3.30,4.24,5.18,6.13,7.07,8.01,8.95};
	double y[]= {5.49,4.08,3.54,2.61,2.09,1.91,1.55,1.47,1.45,1.25};
	double e[]= {0.26,0.12,0.27,0.10,0.15,0.11,0.13,0.07,0.15,0.09};
	int n = sizeof(t)/sizeof(t[0]);

	

	FILE* datastream=fopen("data.txt","w");
	fprintf(datastream, "# t[i] \t y[i] \t e[i]\n");
	for (int i = 0; i < n; i++)
	{
		fprintf(datastream, "%g \t %g \t %g \n", t[i], y[i], e[i]);
	}



	struct experimental_data parameters;
	parameters.n=n;
	parameters.y=y;
	parameters.e=e;
	parameters.t=t;

	int vectorDim=3;
	const gsl_multimin_fminimizer_type* fmintype = 
	gsl_multimin_fminimizer_nmsimplex2;
	
	gsl_multimin_fminimizer* wrkspc = 
	gsl_multimin_fminimizer_alloc(fmintype, vectorDim);

	gsl_vector* startpoint = gsl_vector_alloc(vectorDim);
	gsl_vector_set(startpoint, 0, -5);
	gsl_vector_set(startpoint, 1, 5);
	gsl_vector_set(startpoint, 2, 10);

	gsl_vector* initStepSize = gsl_vector_alloc(vectorDim);
	gsl_vector_set_all(initStepSize, 1.0);

	gsl_multimin_function leastsquares;
	leastsquares.n = vectorDim;
	leastsquares.f = min_fun;
	leastsquares.params=(void*)&parameters;

	gsl_multimin_fminimizer_set(wrkspc,&leastsquares,startpoint,initStepSize);

	int iter=0, flag;

	do
	{
		iter++;
		gsl_multimin_fminimizer_iterate(wrkspc);
		flag = 
		gsl_multimin_test_size(gsl_multimin_fminimizer_size(wrkspc),1e-6);

	}
	while(flag==GSL_CONTINUE && iter<10000);


	double A=gsl_vector_get(wrkspc->x,0);
	double T=gsl_vector_get(wrkspc->x,1);
	double B=gsl_vector_get(wrkspc->x,2);

	printf("Fitting function:\n");
	printf("f(t)=%g*exp(-t/%g)+%g\n",A,T,B);

	fprintf(datastream, "\n\n");
	fprintf(datastream, "# t \t A*exp(-t/T)+B\n");
	for (double t = 0; t < 10; t=t+0.1)
	{
		fprintf(datastream, "%g \t %g\n", t, A*exp(-t/T)+B);
	}

	gsl_vector_free(initStepSize);
	gsl_vector_free(startpoint);
	gsl_multimin_fminimizer_free(wrkspc);
	fclose(datastream);
	return 0;
}