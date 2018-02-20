#ifndef KOMPLEX_H

/*struct komplexVal {double re; double im;}*/
typedef struct {double re; double im;} komplex;

void	komplex_print(char* s, komplex z); /*Given string s and komplex number z, prints string followed by z*/
void	komplex_set(komplex* z, double x, double y); /*Redefine z s.t. z=x+y*i*/
komplex komplex_new(double x, double y); /*Returns x+y*i*/
komplex komplex_add(komplex a, komplex b); /*Returns a+b where a,b are complex numbers*/
komplex komplex_sub(komplex a, komplex b); /*Returns a-b where a,b are complex numbers*/

int		komplex_equal(komplex a, komplex b); /*Returns 1 if equal, 0 otherwise*/
komplex komplex_mul(komplex a, komplex b); /* Returns a*b */
komplex komplex_div(komplex a, komplex b); /* Returns a/b */
komplex komplex_conjugate(komplex z); /*Returns the complex conjugate of z*/
komplex komplex_abs(komplex z); /*Returns absolute value of z*/
komplex komplex_exp(komplex z); /*Complex exponential*/
komplex komplex_cos(komplex z); /*Complex cosine*/
komplex komplex_sin(komplex z); /*Complex sine*/
komplex komplex_sqrt(komplex z); /*Complex principal square root*/

#define KOMPLEX_H
#endif