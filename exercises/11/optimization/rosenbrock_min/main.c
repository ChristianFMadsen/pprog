#include <stdio.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>


double rosenbrock_fun(const gsl_vector* v, void* params)
{
	double x,y;

	x=gsl_vector_get(v,0);
	y=gsl_vector_get(v,1);

	return (1-x)*(1-x)+100*(y-x*x)*(y-x*x);
}

int main()
{
	int vectorDim=2;

	const gsl_multimin_fminimizer_type* fmintype = 
	gsl_multimin_fminimizer_nmsimplex2;
	
	gsl_multimin_fminimizer* wrkspc = 
	gsl_multimin_fminimizer_alloc(fmintype, vectorDim);


	gsl_vector* startpoint = gsl_vector_alloc(vectorDim);
	gsl_vector_set(startpoint, 0, -5);
	gsl_vector_set(startpoint, 1, 5);


	gsl_vector* initStepSize = gsl_vector_alloc(vectorDim);
	gsl_vector_set_all(initStepSize, 1.0);

	gsl_multimin_function rosbro_fun;
	rosbro_fun.n=vectorDim;
	rosbro_fun.f = rosenbrock_fun;
	rosbro_fun.params=NULL;

	gsl_multimin_fminimizer_set(wrkspc,&rosbro_fun,startpoint,initStepSize);

	int iter=0, flag;
	double before_iter, after_iter, deltaf, xFunEst, yFunEst;

	do
	{
		iter++;
		before_iter=gsl_multimin_fminimizer_minimum(wrkspc);
		xFunEst=gsl_vector_get(wrkspc->x,0);
		yFunEst=gsl_vector_get(wrkspc->x,1);
		gsl_multimin_fminimizer_iterate(wrkspc);
		after_iter=gsl_multimin_fminimizer_minimum(wrkspc);
		deltaf=after_iter-before_iter;

		if (deltaf != 0)
		{
			printf("Function estimate at (x,y)=(%g,%g)\n", xFunEst, yFunEst);
		}

		flag = 
		gsl_multimin_test_size(gsl_multimin_fminimizer_size(wrkspc),1e-6);

	} 
	while(flag==GSL_CONTINUE && iter<10000);

	printf("Minimum found at (x,y)=(%g,%g)\n", gsl_vector_get(wrkspc->x,0),gsl_vector_get(wrkspc->x,1));
	printf("In %i iterations\n", iter);

	gsl_vector_free(initStepSize);
	gsl_vector_free(startpoint);
	gsl_multimin_fminimizer_free(wrkspc);
	return 0;
}