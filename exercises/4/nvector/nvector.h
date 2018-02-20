#ifndef NVECTOR_H

typedef struct {int size; double* data;} nvector;

nvector* nvector_alloc(int n); /* Allocates memory for vector of size n*/
void nvector_free(nvector* v); /* Frees up the memory*/
void nvector_set(nvector* v, int i, double value); /*Sets value of i'th place in vector*/
double nvector_get(nvector* v, int i); /*Returns value of i'th place in vector*/
double nvector_dot_product (nvector* u, nvector* v); /*Returns dot product between u and v*/

// Optional
void nvector_print(char* s, nvector* v); //Prints s followed by the vector v
void nvector_set_zero(nvector* v); //Sets all elements in vector equal to zero
int nvector_equal(nvector* u, nvector* v); //Returns 1 if vectors u and v are equal, 0 otherwise.
void nvector_add(nvector* u, nvector* v); // u_i => u_i+v_i
void nvector_sub(nvector* u, nvector* v); // u_i => u_i-v_i
void nvector_scale(nvector* u, double x); //u_i => x*u_i


#define NVECTOR_H
#endif