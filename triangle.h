#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "vector.h"

typedef struct triangle {
	vector *A;
	vector *B;
	vector *C;
	vector *n;
	double reflect;
} triangle;

triangle *new_triangle(vector *A, vector *B, vector *C, double reflect);
void delete_triangle(triangle *t);

#endif // TRIANGLE_H
