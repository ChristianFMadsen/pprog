#include <math.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


/* Jacobi diagonlaization with cyclic sweeps
evals and evecs stores eigenvalues and eigenvectors respectively. 
Only the upper triangular part of the matrix gets updated i.e. destroyed.*/

//Notation: aij = Element on i'th row and j'th column in A

int jacobi(gsl_matrix* A, gsl_matrix* evecs, gsl_vector* evals){
int change, sweepCount=0, n=A->size1; //size1: number of rows

for(int i=0; i<n; i++){
gsl_vector_set(evals,i,gsl_matrix_get(A,i,i)); } //Diagonal elements of A in evals 


gsl_matrix_set_identity(evecs);

do{

	change=0; sweepCount++; int p,q;

	for(p=0; p<n; p++){ //p designates row in A
		for(q=p+1; q<n; q++){ //q designates column in A. 
			double app=gsl_vector_get(evals,p);
			double aqq=gsl_vector_get(evals,q);
			double apq=gsl_matrix_get(A,p,q);
			double phi=0.5*atan2(2*apq, aqq-app); //Angle to zero apq
			double c=cos(phi), s=sin(phi);
			/* eq. 9 in eigen.pdf */
			double app1=c*c*app-2*s*c*apq+s*s*aqq;
			double aqq1=s*s*app+2*s*c*apq+c*c*aqq;

			if(app1!=app || aqq1!=aqq){
				change=1;
				gsl_vector_set(evals,p,app1);
				gsl_vector_set(evals,q,aqq1);
				gsl_matrix_set(A,p,q,0.0);
				

				for(int i=0; i<p; i++){ 
					double aip=gsl_matrix_get(A,i,p); //Update entries above app
					double aiq=gsl_matrix_get(A,i,q); //Update entries above apq
					gsl_matrix_set(A,i,p,c*aip-s*aiq);
					gsl_matrix_set(A,i,q,c*aiq+s*aip); }

				for(int i=p+1; i<q; i++){
					double api=gsl_matrix_get(A,p,i); //Update entries to the right of app
					double aiq=gsl_matrix_get(A,i,q); //Update entries under apq
					gsl_matrix_set(A,p,i,c*api-s*aiq);
					gsl_matrix_set(A,i,q,c*aiq+s*api); }

				for(int i=q+1; i<n; i++){
					double api=gsl_matrix_get(A,p,i); //Update entries to the right of apq
					double aqi=gsl_matrix_get(A,q,i); //Update entries to the right of aqq
					gsl_matrix_set(A,p,i,c*api-s*aqi);
					gsl_matrix_set(A,q,i,c*aqi+s*api); }

				for(int i=0; i<n; i++){ //Fill eigenvectors into evecs
					double evecs_ip=gsl_matrix_get(evecs, i, p);
					double evecs_iq=gsl_matrix_get(evecs, i, q);
					gsl_matrix_set(evecs, i, p, c*evecs_ip-s*evecs_iq);
					gsl_matrix_set(evecs, i, q, c*evecs_iq+s*evecs_ip); }	
			}
		}
	}
}while(change!=0);


return sweepCount; 
}

int jacobi_eig_ascending(gsl_matrix* A, gsl_vector* evalsVec, int k){

int n=A->size1;
assert(k <= n);	
assert(evalsVec->size == k);
int change, p, q;
gsl_vector* evals = gsl_vector_calloc(n);

for(int i=0; i<n; i++){
gsl_vector_set(evals,i,gsl_matrix_get(A,i,i)); }


for(p=0; p<k; p++){

do{
	change=0;
		for(q=p+1; q<n; q++){ //q designates column in A. 
			double app=gsl_vector_get(evals,p);
			double aqq=gsl_vector_get(evals,q);
			double apq=gsl_matrix_get(A,p,q);
			double phi=0.5*atan2(2*apq, aqq-app); 
			double c=cos(phi), s=sin(phi);
			double app1=c*c*app-2*s*c*apq+s*s*aqq;
			double aqq1=s*s*app+2*s*c*apq+c*c*aqq;

			if(app1!=app || aqq1!=aqq){
				change=1;
				gsl_vector_set(evals,p,app1);
				gsl_vector_set(evals,q,aqq1);
				gsl_matrix_set(A,p,q,0.0);
				
				for(int i=p+1; i<q; i++){
					double api=gsl_matrix_get(A,p,i); //Update entries to the right of app
					double aiq=gsl_matrix_get(A,i,q); //Update entries under apq
					gsl_matrix_set(A,p,i,c*api-s*aiq);
					gsl_matrix_set(A,i,q,c*aiq+s*api); }

				for(int i=q+1; i<n; i++){
					double api=gsl_matrix_get(A,p,i); //Update entries to the right of apq
					double aqi=gsl_matrix_get(A,q,i); //Update entries to the right of aqq
					gsl_matrix_set(A,p,i,c*api-s*aqi);
					gsl_matrix_set(A,q,i,c*aqi+s*api); }
			}
		}
	}while(change!=0); 
}

for(int i=0; i<k; i++){
	gsl_vector_set(evalsVec, i, gsl_vector_get(evals,i)); }


gsl_vector_free(evals);
return 0;
}


int jacobi_eig_descending(gsl_matrix* A, gsl_vector* evalsVec, int k){

int n=A->size1;
assert(k <= n);	
assert(evalsVec->size == k);
int change, p, q;

gsl_vector* evals = gsl_vector_calloc(n);
for(int i=0; i<n; i++){
gsl_vector_set(evals,i,gsl_matrix_get(A,i,i)); }


for(p=0; p<k; p++){

do{
	change=0;
		for(q=p+1; q<n; q++){ //q designates column in A. 
			double app=gsl_vector_get(evals,p);
			double aqq=gsl_vector_get(evals,q);
			double apq=gsl_matrix_get(A,p,q);
			double phi=0.5*atan2(2*apq, -aqq+app); //Change sign in denominator
			double c=cos(phi), s=-sin(phi); //Change rotation matrix, this is the transpose of the original rotation matrix.
			double app1=c*c*app-2*s*c*apq+s*s*aqq;
			double aqq1=s*s*app+2*s*c*apq+c*c*aqq;

			if(app1!=app || aqq1!=aqq){
				change=1;
				gsl_vector_set(evals,p,app1);
				gsl_vector_set(evals,q,aqq1);
				gsl_matrix_set(A,p,q,0.0);
				
				for(int i=p+1; i<q; i++){
					double api=gsl_matrix_get(A,p,i); //Update entries to the right of app
					double aiq=gsl_matrix_get(A,i,q); //Update entries under apq
					gsl_matrix_set(A,p,i,c*api-s*aiq);
					gsl_matrix_set(A,i,q,c*aiq+s*api); }

				for(int i=q+1; i<n; i++){
					double api=gsl_matrix_get(A,p,i); //Update entries to the right of apq
					double aqi=gsl_matrix_get(A,q,i); //Update entries to the right of aqq
					gsl_matrix_set(A,p,i,c*api-s*aqi);
					gsl_matrix_set(A,q,i,c*aqi+s*api); }
			}
		}
	}while(change!=0); 
}

for(int i=0; i<k; i++){
	gsl_vector_set(evalsVec, i, gsl_vector_get(evals,i)); }


gsl_vector_free(evals);
return 0;
}