#include "triangle.h"
#include "vector.h"
#include <stdlib.h>

triangle *new_triangle(vector *A, vector *B, vector *C, double reflect) {
	triangle *new = (triangle *) calloc(1, sizeof(triangle));
	new->A = A;
	new->B = B;
	new->C = C;
	new->reflect = reflect;

	new->n = new_zvector();
	vector *ab = new_zvector();
	vector *ac = new_zvector();
	sub(B, A, ab);
	sub(C, A, ac);
	cross(ab, ac, new->n);
	free(ab);
	free(ac);

	return new;
}

void delete_triangle(triangle *t) {
	free(t->A);
	free(t->B);
	free(t->C);
	free(t->n);

	free(t);
}
