#include "vector.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

vector *new_vector(double x, double y, double z) {
	vector *new = (vector *) calloc(1, sizeof(vector));
	new->x = x;
	new->y = y;
	new->z = z;

	return new;
}

vector *new_zvector() {
	return (vector *) calloc(1, sizeof(vector));
}

vector *cpy_vector(const vector *v) {
	vector *new = (vector *) calloc(1, sizeof(vector));
	new->x = v->x;
	new->y = v->y;
	new->z = v->z;

	return new;
}

void set(vector *v, double x, double y, double z) {
	v->x = x;
	v->y = y;
	v->z = z;
}

void assign(vector *v, const vector *u) {
	v->x = u->x;
	v->y = u->y;
	v->z = u->z;
}

void add(const vector *l, const vector *r, vector *res) {
	set(res, l->x + r->x,
	         l->y + r->y,
		 l->z + r->z);
}

void sub(const vector *l, const vector *r, vector *res) {
	set(res, l->x - r->x,
	         l->y - r->y,
		 l->z - r->z);
}

void mul(const vector *l, double r, vector *res) {
	set(res, l->x * r,
	         l->y * r,
		 l->z * r);
}

void move(vector *v, const vector *oth) {
	v->x += oth->x;
	v->y += oth->y;
	v->z += oth->z;
}

void scale(vector *v, double factor) {
	v->x *= factor;
	v->y *= factor;
	v->z *= factor;
}

double dot(const vector *l, const vector *r) {
	return l->x * r->x +
	       l->y * r->y +
	       l->z * r->z;
}

void cross(const vector *l, const vector *r, vector *res) {
	res->x = l->y * r->z - l->z * r->y;
	res->y = l->z * r->x - l->x * r->z;
	res->z = l->x * r->y - l->y * r->x;
}

double mixed(const vector *a, const vector *b, const vector *c) {
	return a->x * (b->y * c->z - b->z * c->y) -
	       a->y * (b->x * c->z - b->z * c->x) +
	       a->z * (b->x * c->y - b->y * c->x);
}

double len_sq(const vector *v) {
	return v->x * v->x +
	       v->y * v->y +
	       v->z * v->z;
}

double len(const vector *v) {
	return sqrt(len_sq(v));
}

void dump(const vector *v) {
	printf("(%g, %g, %g)\n", v->x, v->y, v->z);
}

