#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "cell.h"

#define MAX_ITERATIONS 10000

// ----------------- UTILS -----------------

double rand_double(double high) {
    return ((double)rand() * high) / (double)RAND_MAX;
}

double distance(const Point* p1, const Point* p2) {
    return sqrt(pow(p1->x - p2->x, 2) + pow(p1->y - p2->y, 2));
}

double distance_pbc(const Point* p1, const Point* p2, const GenParams *params) {
    double dx = fabs(p1->x - p2->x);
    double dy = fabs(p1->y - p2->y);
    double dist_x = fmin(dx, params->Lx - dx);
    double dist_y = fmin(dy, params->Ly - dy);
    return sqrt(pow(dist_x, 2) + pow(dist_y, 2));
}

int random_gen_point(const GenParams * params, Point* points, size_t t_idx, bool pbc) {
    for (size_t iter_idx = 0; iter_idx < MAX_ITERATIONS; iter_idx++) {
        points[t_idx].x = rand_double(params->Lx);
        points[t_idx].y = rand_double(params->Ly);
        size_t p_idx;
        for (p_idx = 0; p_idx < t_idx; p_idx++) {
            double dist = pbc ? distance_pbc(points + p_idx, points + t_idx, params) : distance(points + p_idx, points + t_idx);
            if (dist < params->sigma) break;
        }
        if (p_idx == t_idx) return 0;
    }
    return 1;
}

// ----------------- MAIN FUNCTIONS -----------------

int write_pdb(const char* filename, const Point* points, size_t N) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) return 1;
    for (size_t idx = 0; idx < N; idx++)
        fprintf(file, "ATOM  %5ld  N   NONE    1    %8.3f%8.3f%8.3f  1.00  0.00\n", idx + 1, points[idx].x, points[idx].y, 0.0);
    fclose(file);
    return 0;
}

int random_gen(const char* filename, const GenParams *params, bool pbc) {
    Point *points = (Point*)malloc(params->N * sizeof(Point));
    size_t num_generated = params->N;
    for (size_t idx = 0; idx < params->N; idx++) {
        if(random_gen_point(params, points, idx, pbc)) {
            num_generated = idx;
            break;
        }
    }
    if (num_generated != params->N) {
        printf("Could not generate %d points\n", params->N);
    }
    int write_code = write_pdb(filename, points, num_generated);
    free(points);
    return write_code;
}

int square_gen(const char* filename, const GenParams *params, bool pbc) {
    unsigned int N_x = (unsigned int) (params->Lx / params->sigma) + 1;
    unsigned int N_y = (unsigned int) (params->Ly / params->sigma) + 1;
    unsigned int N = N_x * N_y;
    Point *points = (Point*)malloc(N * sizeof(Point));
    for (unsigned int x_idx = 0; x_idx < N_x; x_idx++) {
        for (unsigned int y_idx = 0; y_idx < N_y; y_idx++) {
            points[x_idx * N_y + y_idx].x = x_idx * params->sigma;
            points[x_idx * N_y + y_idx].y = y_idx * params->sigma;
        }
    }
    int write_code = write_pdb(filename, points, N);
    free(points);
    return write_code;
}

int hexagonal_gen(const char* filename, const GenParams *params, bool pbc) {
    unsigned int N_x = (unsigned int) (params->Lx / params->sigma / (sqrt(3)/2));
    unsigned int N_y = (unsigned int) (params->Ly / params->sigma);
    unsigned int N = N_x * N_y;
    Point *points = (Point*)malloc(N * sizeof(Point));
    for (unsigned int x_idx = 0; x_idx < N_x; x_idx++) {
        for (unsigned int y_idx = 0; y_idx < N_y; y_idx++) {
            points[x_idx * N_y + y_idx].x = x_idx * params->sigma * sqrt(3) / 2;
            points[x_idx * N_y + y_idx].y = y_idx * params->sigma;
            if (x_idx % 2 == 1) points[x_idx * N_y + y_idx].y += params->sigma / 2;
        }
    }
    int write_code = write_pdb(filename, points, N);
    free(points);
    return write_code;
}
