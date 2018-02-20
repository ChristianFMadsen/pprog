#include "nvector.h"
#include <stdio.h>
#include <stdlib.h>

nvector* nvector_alloc(int n){
	nvector* v = malloc(sizeof(nvector)); /*wut*/
	/* printf("%lu\n", malloc(sizeof(nvector))); */
	(*v).size = n;
	(*v).data = malloc(n*sizeof(double));
	if (v==NULL){
		fprintf(stderr, "Error in nvector_alloc\n");
	}
	return v;
}

void nvector_free(nvector* v){
	free((*v).data); /* redundant?*/
	free(v);
	v=NULL;
}

void nvector_set(nvector* v, int i, double value){
	(*v).data[i]=value;
}

double nvector_get(nvector* v, int i){
	return (*v).data[i];
}

double nvector_dot_product (nvector* u, nvector* v){
	double result=0.0;
	if((*u).size == (*v).size){
		for(int i=0; i<(*u).size; i++){
			result=result+(*u).data[i]*(*v).data[i];
		}
		printf("%g\n", result);
	}
	else {
		printf("Dimensions does not match\n");
	}
	return result;
}

void nvector_print(char* s, nvector* v){
	printf("%s\n", s);
	int i=0;
	printf("("); //Left parenthesis
	while (i<(*v).size-2){
		printf("%g, ",(*v).data[i]); //Prints vector contents separated by ,
		i++;
	}
	printf("%g)\n", (*v).data[i+1]); //Prints last entry in vector and the right parenthesis 
}

void nvector_set_zero(nvector* v){
	for (int i=0; i<(*v).size; i++){
		(*v).data[i]=0;
	}
}

int nvector_equal(nvector* u, nvector* v){
	if((*u).size == (*v).size){ //Check if same dimensions
		for (int i=0; i<(*u).size; i++){ //Compare contents of vectors
			if((*u).data[i]==(*v).data[i]){ //If i'th element of vectors equal do nothing
			}
			else{ //If not equal return 0 and break
				return 0;
			}
		}
		return 1; //Return 1 since all elements were equal
	}
	else{ //If not same dimensions return 0.
		return 0;
	}
}

void nvector_add(nvector* u, nvector* v){
	if((*u).size == (*v).size){
		for(int i=0; i<(*u).size; i++){
			(*u).data[i]=(*u).data[i]+(*v).data[i];
		}
	}
	else{
		printf("Dimensions does not match\n");
	}
}

void nvector_sub(nvector* u, nvector* v){
	if((*u).size == (*v).size){
		for(int i=0; i<(*u).size; i++){
			(*u).data[i]=(*u).data[i]-(*v).data[i];
		}
	}
	else{
		printf("Dimensions does not match\n");
	}	
}

void nvector_scale(nvector* u, double x){
	for(int i=0; i<(*u).size; i++){
		(*u).data[i]=(*u).data[i]*x;
	}
}