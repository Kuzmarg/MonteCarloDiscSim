#ifndef UTILS_H
#define UTILS_H

#include <stdbool.h>
#include <curand_kernel.h>
#include "types.h"

////////////////////////////////////////////////////////////////////////////////

__host__ __device__ void quaternion_to_angle(double qz, double qw, double *sin_t, double *cos_t);

__host__ __device__ void rotate_point(double x, double y, double cos_t, double sin_t, double *x_rot, double *y_rot);

__host__ __device__ void get_square_vertices(double x, double y, double s, double cos_t, double sin_t, double *vertices);

////////////////////////////////////////////////////////////////////////////////

__host__ double rand_double(double high);

__host__ int rand_int(int high);

////////////////////////////////////////////////////////////////////////////////

__global__ void rand_init_kernel(curandState *state);

__device__ double rand_double_cuda(double high, curandState *state);

__device__ int rand_int_cuda(int high, curandState *state);

////////////////////////////////////////////////////////////////////////////////

__host__ __device__ double distance(const Particle* p1, const Particle* p2, const Config* config);

__host__ __device__ double distance_patch(const Patch* p1, const Patch* p2, const Config* config);

////////////////////////////////////////////////////////////////////////////////

__host__ int write_xyz(const char* filename, const Config* config, const CellLinkedGrid* cll);

#endif // UTILS_H
