#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include "types.h"

void quaternion_to_angle(double qz, double qw, double *sin_t, double *cos_t);

void rotate_point(double x, double y, double cos_t, double sin_t, double *x_rot, double *y_rot);

void get_square_vertices(double x, double y, double s, double cos_t, double sin_t, double *vertices);

double rand_double(double high);

double distance(const Particle* p1, const Particle* p2, const Grid* grid);

int write_pdb(const char* filename, const Grid* grid);

int write_xyz(const char* filename, const Grid* grid);

#endif // UTILS_H
