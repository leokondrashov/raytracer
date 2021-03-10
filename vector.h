#ifndef VECTOR_H
#define VECTOR_H

typedef struct vector {
	double x;
	double y;
	double z;
} vector;

vector *new_vector(double x, double y, double z);
vector *new_zvector();
vector *cpy_vector(const vector *v);
void set(vector *v, double x, double y, double z);
void assign(vector *v, const vector *u);

void add(const vector *l, const vector *r, vector *res);
void sub(const vector *l, const vector *r, vector *res);
void mul(const vector *l, double r, vector *res);

double dot(const vector *l, const vector *r);
void cross(const vector *l, const vector *r, vector *res);
double mixed(const vector *a, const vector *b, const vector *c);

double len_sq(const vector *v);
double len(const vector *v);

void dump(const vector *v);

#endif // VECTOR_H
