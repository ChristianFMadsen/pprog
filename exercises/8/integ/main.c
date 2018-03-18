#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

//Functions to be used in gsl_function calls
double log_integrand (double x, void* params) {
  double f = log(x) / sqrt(x);
  return f;
}

double norm (double x, void* params) {
	double alpha = *(double*)params;
	double f = exp(-alpha*x*x);
	return f;
}

double hamiltonian (double x, void* params) {
	double alpha = *(double*)params;
	double f = (-alpha*alpha*x*x/2.0+alpha/2.0+x*x/2.0)*exp(-alpha*x*x);
	return f;
}

int main(){
	//Part 1:
	double expected_result = -4.0; //Obtained from WolframAlpha
	printf("Exact result of log(x)/sqrt(x) integral: %.18f\n", expected_result);

	double result, error, a, b, epsabs, epsrel, limit;
	a=0; //Integrate from
	b=1; //Integrate to
	epsabs=1e-7;
	epsrel=1e-7;
	limit = 1000.0;

	gsl_integration_workspace* w1 = gsl_integration_workspace_alloc(limit); //create workspace
	gsl_function logF; //create function
	logF.function=&log_integrand;
	logF.params=NULL;

	gsl_integration_qags(&logF,a,b,epsabs,epsrel,limit,w1,&result, &error); //Integrate
	printf("result = %.18f\n", result);
	printf("estimated error = %.18f\n", error);
	printf("actual error = %.18f\n", result-expected_result);

	//Part 2
	double alpha, result_norm, error_norm, result_hamiltonian, error_hamiltonian;

	gsl_function Norm;
	Norm.function=&norm;
	Norm.params=&alpha; 

	gsl_function Hamiltonian;
	Hamiltonian.function=&hamiltonian;
	Hamiltonian.params=&alpha;

	FILE* datastream = fopen("energies.txt","w");

	for (alpha=0.1; alpha<10; alpha=alpha+0.1){
		gsl_integration_qagi(&Norm,epsabs,epsrel,limit,w1,&result_norm,&error_norm);		
		gsl_integration_qagi(&Hamiltonian,epsabs,epsrel,limit,w1,&result_hamiltonian,&error_hamiltonian);
		fprintf(datastream, "%.1f \t %.18f\n", alpha, result_hamiltonian/result_norm);
	}



	fclose(datastream);
	gsl_integration_workspace_free(w1);
	return 0;
}