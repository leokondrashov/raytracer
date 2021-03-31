#include "vector.h"
#include "bmp.h"
#include "triangle.h"
#include "textProcessor/textProcessor.h"
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define WIDTH 1024
#define HEIGHT 1024

int collision(const vector *ray, const vector *origin, const triangle *t, vector *P, double *r) {
	static vector ABC[3];
	sub(t->A, origin, &(ABC[0]));
	sub(t->B, origin, &(ABC[1]));
	sub(t->C, origin, &(ABC[2]));
	double ab = mixed(&(ABC[0]), &(ABC[1]), ray);
	double bc = mixed(&(ABC[1]), &(ABC[2]), ray);
	double ca = mixed(&(ABC[2]), &(ABC[0]), ray);

	int collided = //(ab >= 0 && bc >= 0 && ca >= 0) ||
	               (ab <= 0 && bc <= 0 && ca <= 0);
	double proj = dot(&(ABC[0]), t->n) / dot(ray, t->n);
	if (!collided || proj < 0)
		return 0;

	*r = len_sq(ray) * proj * proj;

	if (P != NULL) {
		assign(P, ray);
		scale(P, proj);
		move(P, origin);
	}

	return 1;
}

double pixel_color(const vector *ray, const vector *origin, triangle *const *t, 
                   int t_n, vector *const *l, int l_n) {
	double min_dist = -1;
	int min_i = -1;
	static vector min_P;

	static vector P;
	double r = 0;
	int collided = 0;
	int skip = -1;

	double weight = 1;
	double color = 0;
	double cur_light = 0;

	static vector new_ray;
	static vector o;
	assign(&o, origin);
	static vector R;
	assign(&R, ray);

	do {
		collided = 0;
		min_dist = -1;
		min_i = -1;
		r = 0;

		for (int i = 0; i < t_n; ++i) {
			if (i == skip) 
				continue;

			if (collision(&R, &o, t[i], &P, &r) && 
			    (r < min_dist || min_dist == -1)) {
				collided = 1;
				min_dist = r;
				min_i = i;
				assign(&min_P, &P);
			}
		}

		if (collided) {
			cur_light = 0;

			for (int i = 0; i < l_n; ++i) {
				sub(&min_P, l[i], &new_ray);
				int min_j = -1;
				min_dist = -1;
				for (int j = 0; j < t_n; ++j) {
					if (collision(&new_ray, l[i], t[j], NULL, &r) &&
					    (r < min_dist || min_dist == -1)) {
						min_dist = r;
						min_j = j;
					}
				}

				double cos = -dot(t[min_i]->n, &new_ray) /
					     len(t[min_i]->n) / len(&new_ray);
				cur_light += (cos > 0 && min_j == min_i) ? cos : 0;
			}

			assign(&new_ray, t[min_i]->n);
			scale(&new_ray, -2 * dot(&R, &new_ray) / len_sq(&new_ray));
			move(&new_ray, &R);
			color += cur_light * (1 - t[min_i]->reflect) * weight;
			weight *= t[min_i]->reflect;

			assign(&R, &new_ray);
			assign(&o, &min_P);
			skip = min_i;
		}
	} while (collided);

	return color;
}

void parse(const char *file, triangle ***t, vector ***l, int *t_n, int *l_n) {
	char **txt = readTextFromFile(file);
	int i = 0;

	sscanf(txt[i++], "%d", t_n);

	*t = (triangle **) calloc(*t_n, sizeof(triangle *));
	double ax = 0, ay = 0, az = 0;
	double bx = 0, by = 0, bz = 0;
	double cx = 0, cy = 0, cz = 0;
	double reflect = 0;
	for (int j = 0; j < *t_n; ++j) {
		sscanf(txt[i++], "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &ax, &ay, &az,
		                                                            &bx, &by, &bz,
		                                                            &cx, &cy, &cz,
		                                                            &reflect);
		(*t)[j] = new_triangle(new_vector(ax, ay, az),
		                       new_vector(bx, by, bz), 
		                       new_vector(cx, cy, cz), reflect);
	}

	sscanf(txt[i++], "%d", l_n);

	*l = (vector **) calloc(*l_n, sizeof(vector *));
	double x = 0, y = 0, z = 0;
	for (int j = 0; j < *l_n; ++j) {
		sscanf(txt[i++], "%lg %lg %lg", &x, &y, &z);
		(*l)[j] = new_vector(x, y, z);
	}
	
	free(*txt);
	free(txt);
}

int main(int argc, char **argv) {
//	MPI_Init(&argc, &argv);

	int r = 0, s = 1;
//	MPI_Comm_rank(MPI_COMM_WORLD, &r);
//	MPI_Comm_size(MPI_COMM_WORLD, &s);
	int n = HEIGHT / s;
	int is = r * n, ie = (r + 1) * n;
	if (r == s - 1)
		ie = HEIGHT;

	triangle **t = NULL;
	int t_n = 0;
	vector **l = NULL;
	int l_n = 0;
	parse("triangles.txt", &t, &l, &t_n, &l_n);

	char *data = (char *) calloc(WIDTH * HEIGHT, sizeof(char));
	vector *pixel = new_zvector();
	vector *ray = new_zvector();
	vector *origin = new_vector(0, 0, 0);
	vector *P = new_zvector();
	for (int i = is; i < ie; ++i) {
		for (int j = 0; j < WIDTH; ++j) {
			set(pixel, -1 + 2.0 / WIDTH * j, -1 + 2.0 / HEIGHT * i, 1);
			sub(pixel, origin, ray);
			data[i * WIDTH + j] = 255 * pixel_color(ray, origin, t, t_n, l, l_n) / l_n;
		}
	}

	for (int i = 0; i < t_n; ++i)
		delete_triangle(t[i]);
	free(t);

	free(pixel);
	free(ray);
	free(origin);
	free(P);

	for (int i = 0; i < l_n; ++i)
		free(l[i]);
	free(l);

	if (r == s - 1) {
//		for (int i = 0; i < s - 1; ++i)
//			MPI_Recv(data + i * n * WIDTH, n * WIDTH, MPI_CHAR, i, 0, 
//			         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	} else {
//		MPI_Send(data + is * WIDTH, (ie - is) * WIDTH, MPI_CHAR, s - 1, 0, 
//		         MPI_COMM_WORLD);
	}

	if (r == s - 1)
		dump_bmp("1.bmp", data, HEIGHT, WIDTH);

	free(data);
//	MPI_Finalize();
}
