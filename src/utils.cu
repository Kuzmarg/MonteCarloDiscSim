#include "utils.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <curand_kernel.h>

__host__ __device__ void get_square_vertices(double x, double y, double s, double cos_t, double sin_t, double *vertices) {
    double x_rot, y_rot;
    for (int i = 0; i < 4; i++) {
        rotate_point((i % 2) ? -s/2 : +s/2, (i < 2) ? s/2 : -s/2, cos_t, sin_t, &x_rot, &y_rot);
        vertices[2*i] = x + x_rot;
        vertices[2*i + 1] = y + y_rot;
    } 
}

__host__ __device__ void quaternion_to_angle(double qz, double qw, double *sin_t, double *cos_t) {
    *cos_t = 2*qw*qw - 1;
    *sin_t = 2*qw*qz;
}

__host__ __device__ void rotate_point(double x, double y, double cos_t, double sin_t, double *x_rot, double *y_rot) {
    *x_rot = x*cos_t - y*sin_t;
    *y_rot = x*sin_t + y*cos_t;
}

////////////////////////////////////////////////////////////////////////////////

__host__ double rand_double(double high) {
    return ((double)rand() * high) / (double)RAND_MAX;
}

__host__ int rand_int(int high) {
    return rand() % high;
}

////////////////////////////////////////////////////////////////////////////////

__global__ void rand_init_kernel(curandState *state) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    curand_init(1234, idx, 0, &state[idx]);
}

__device__ double rand_double_cuda(double high, curandState *state) {
    return curand_uniform(state) * high;
}

__device__ int rand_int_cuda(int high, curandState *state) {
    return (int)(curand_uniform(state) * high);
}

////////////////////////////////////////////////////////////////////////////////

__host__ __device__ double distance(const Particle* p1, const Particle* p2, const Config* config) {
    double dx = fmin(fabs(p1->x - p2->x), fabs(config->Lx - fabs(p1->x - p2->x)));
    double dy = fmin(fabs(p1->y - p2->y), fabs(config->Ly - fabs(p1->y - p2->y)));
    return sqrt(pow(dx, 2) + pow(dy, 2));
}

__host__ __device__ double distance_patch(const Patch *p1, const Patch *p2, const Config *config) {
    double dx = fmin(fabs(p1->x - p2->x), fabs(config->Lx - fabs(p1->x - p2->x)));
    double dy = fmin(fabs(p1->y - p2->y), fabs(config->Ly - fabs(p1->y - p2->y)));
    return sqrt(pow(dx, 2) + pow(dy, 2));
}

////////////////////////////////////////////////////////////////////////////////

__host__ int write_xyz(const char* filename, const Config* config, const CellLinkedGrid* cll) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) return 1;
    // 92 is the size of a line for each particle, 100 is for the header
    char *file_string = (char *)malloc((92*config->N*(config->num_patches + 1) + 100) * sizeof(char));
    size_t offset = 0;

    char type = config->type == CIRCLE ? 'C' : 'S';
    offset += sprintf(file_string + offset, "%d\n", config->N + config->N * config->num_patches);
    offset += sprintf(file_string + offset, "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3\n");
    for (size_t idx = 0; idx < config->N; idx++) {
        offset += sprintf(file_string + offset, "%c %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
            type,
            cll->particles[idx].x,
            cll->particles[idx].y,
            0.0,
            0.0,
            0.0,
            cll->particles[idx].qz,
            cll->particles[idx].qw,
            config->size/2,
            config->size/2,
            0.2
        );

        // Write patches
        for (size_t j = 0; j < config->num_patches; j++) {
            double cos_t = 2*cll->particles[idx].qw*cll->particles[idx].qw - 1;
            double sin_t = 2*cll->particles[idx].qw*cll->particles[idx].qz;
            double x_rot, y_rot;
            rotate_point(config->patch_coordinates[j].x, config->patch_coordinates[j].y, cos_t, sin_t, &x_rot, &y_rot);
            offset += sprintf(file_string + offset, "P %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                cll->particles[idx].x + x_rot,
                cll->particles[idx].y + y_rot,
                0.0,
                0.0,
                0.0,
                cll->particles[idx].qz,
                cll->particles[idx].qw,
                config->patch_size/2,
                config->patch_size/2,
                0.2
            );
        }
    }
    fprintf(file, "%s", file_string);
    fclose(file);
    return 0;
}
