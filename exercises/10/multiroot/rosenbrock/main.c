#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

int rosenbrock_fun(const gsl_vector* x, void* params, gsl_vector* f)
{
	const double x0=gsl_vector_get(x,0); //x
	const double x1=gsl_vector_get(x,1); //y

	const double y0=2*x0-2-400*x0*(x1-x0*x0); //f_1(x,y)
	const double y1=200*(x1-x0*x0); //f_2(x,y)

	gsl_vector_set(f,0,y0);
	gsl_vector_set(f,1,y1);

	return GSL_SUCCESS;
}

int main()
{
	int vectorDim = 2;
	const gsl_multiroot_fsolver_type* solvertype = gsl_multiroot_fsolver_hybrids;
	gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc(solvertype,vectorDim);


	gsl_multiroot_function ros_fun;
	ros_fun.f=rosenbrock_fun;
	ros_fun.n=vectorDim;
	ros_fun.params=NULL;

	double startx, starty;
	startx=2;
	starty=0;

	gsl_vector* startpoint = gsl_vector_alloc(vectorDim);
	gsl_vector_set(startpoint,0,startx);
	gsl_vector_set(startpoint,1,starty);
	gsl_multiroot_fsolver_set(solver,&ros_fun,startpoint);

	printf("Start guess (x,y) = (%g,%g)\n", startx,starty);

	double before_iterx,before_itery,after_iterx,after_itery,deltax,deltay;
	int flag,iter;
	iter=0;
	do
	{
		iter++;
		before_iterx = gsl_vector_get(solver->f, 0);
		before_itery = gsl_vector_get(solver->f, 1);
		gsl_multiroot_fsolver_iterate(solver);
		after_iterx = gsl_vector_get(solver->f, 0);
		after_itery = gsl_vector_get(solver->f, 1);
		deltax=after_iterx-before_iterx;
		deltay=after_itery-before_itery;

		if(deltax != 0 || deltay != 0)
		{
			printf("The gradient was computed at (x,y)=(%g,%g)\n", gsl_vector_get(solver->x,0),gsl_vector_get(solver->x,1));
		}

		flag = gsl_multiroot_test_residual(solver->f,1e-12);
		
	}
	while(flag==GSL_CONTINUE);
	double resultx = gsl_vector_get(solver->x,0);
	double resulty = gsl_vector_get(solver->x,1);

	printf("Root found at: (x,y)=(%g,%g)\n",resultx,resulty);
	printf("Root found in %i iterations.\n", iter);



	gsl_vector_free(startpoint);
	gsl_multiroot_fsolver_free(solver);
	return 0;
}