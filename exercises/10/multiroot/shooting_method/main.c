#include<stdio.h>
#include<math.h>
#include <gsl/gsl_vector.h>
#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_errno.h>

double Fepsi(double eps, double r);

int Hfun(const gsl_vector* x, void* params, gsl_vector* f)
{
	double eps=gsl_vector_get(x,0);
	double rmax=8;
	gsl_vector_set(f,0,Fepsi(eps,rmax));
	return GSL_SUCCESS;
}

int main()
{
	int vectorDim=1;
	gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids,vectorDim);

	gsl_multiroot_function FUN;
	FUN.f=Hfun;
	FUN.n=vectorDim;
	FUN.params=NULL;

	gsl_vector* startpoint = gsl_vector_alloc(vectorDim);
	gsl_vector_set(startpoint,0,-1);

	gsl_multiroot_fsolver_set(solver, &FUN, startpoint);

	int flag, iter;
	iter=0;

	do
	{
		iter++;
		gsl_multiroot_fsolver_iterate(solver);
		flag = gsl_multiroot_test_residual(solver->f,1e-10);
	}
	while(flag == GSL_CONTINUE);
	double eps=gsl_vector_get(solver->x,0);
	printf("Calculated eps: %g\n", eps);

	
	FILE* datastream=fopen("data.txt", "w");
	double rmax=8;
	fprintf(datastream,"r \t Fepsi(eps,r) \t exact:\n");

	for(double r=0; r<=rmax; r=r+rmax/100)
	{
		fprintf(datastream,"%.6g \t %.6g \t %.6g\n", r, Fepsi(eps,r), r*exp(-r));
	}
	

	fclose(datastream);
	gsl_multiroot_fsolver_free(solver);
	gsl_vector_free(startpoint);
	return 0;
}