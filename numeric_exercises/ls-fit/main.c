#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<stdio.h>
#include "qr-ls.h"

int main(){
//Data:
double x[]={0.1, 1.33, 2.55, 3.78, 5.0, 6.22, 7.45, 8.68, 9.9};
double y[]={-15.3, 0.32, 2.45, 2.75, 2.27, 1.35, 0.157, -1.23, -2.75};
double errors[]={1.04, 0.594, 0.983, 0.998, 1.11, 0.398, 0.535, 0.968, 0.478};
double x1[] = {0.100,0.145,0.211,0.307,0.447,0.649,0.944,1.372,1.995,2.900};
double y1[] = {12.644,9.235,7.377,6.460,5.555,5.896,5.673,6.964,8.896,11.355};
double errors1[] = {0.858,0.359,0.505,0.403,0.683,0.605,0.856,0.351,1.083,1.002};

int n=sizeof(x)/sizeof(x[0]); //x, y, errors vector dimension
int colDim=3;


// Vectors/matrices for leastsquares function
gsl_vector* xVec=gsl_vector_calloc(n);
gsl_vector* yVec=gsl_vector_calloc(n);
gsl_vector* errorsVec=gsl_vector_calloc(n);

for(int i=0; i<n; i++){
	gsl_vector_set(xVec, i, x[i]);
	gsl_vector_set(yVec, i, y[i]);
	gsl_vector_set(errorsVec, i, errors[i]); 
}
gsl_matrix* S = gsl_matrix_alloc(colDim, colDim);
gsl_vector* c = gsl_vector_alloc(colDim);
leastsquares(funs, colDim, xVec, yVec, errorsVec, c, S);

gsl_vector* coeffErrors = gsl_vector_calloc(colDim); 
for(int k=0; k<colDim; k++){
	double Skk=gsl_matrix_get(S,k,k);
	gsl_vector_set(coeffErrors,k,sqrt(Skk)); }

//Function to calculate fit according to eq. 5
double fit(double x){
	double res=0;
	for(int k=0;k<colDim;k++){
		res=res+gsl_vector_get(c,k)*funs(k,x); }
	return res;
	}
	
//Function to calculate the uncertainty
double uncfit(double x, gsl_matrix* S){
	double res=0;
	for(int i=0; i<colDim; i++){
		for(int j=0; j<colDim; j++){
			res=res+funs(i,x)*funs(j,x)*gsl_matrix_get(S,i,j); }
	}
return sqrt(res);
}

//index 0: data
printf("#x \t y \t errors\n");
for(int i=0;i<n;i++){ 
	printf("%g \t %g \t %g\n",x[i],y[i],errors[i]); }
printf("\n\n");

double intervals=200;
double deltaX=(x[n-1]-x[0])/intervals;
double j=x[0];

//index 1: fit data
printf("#x \t fit \t fit+uncertainties \t fit-uncertainties\n");
while(j<x[n-1]){
	printf("%g \t %g \t %g \t %g \n",j,fit(j),fit(j)+uncfit(j,S),fit(j)-uncfit(j,S));
	j=j+deltaX; }


gsl_vector_free(xVec); gsl_vector_free(yVec); gsl_vector_free(errorsVec); 
gsl_matrix_free(S); gsl_vector_free(c); gsl_vector_free(coeffErrors);
return 0;
}