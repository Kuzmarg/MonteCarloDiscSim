#include "utils.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void get_square_vertices(double x, double y, double s, double cos_t, double sin_t, double *vertices) {
    double x_rot, y_rot;
    for (int i = 0; i < 4; i++) {
        rotate_point((i % 2) ? -s/2 : +s/2, (i < 2) ? s/2 : -s/2, cos_t, sin_t, &x_rot, &y_rot);
        vertices[2*i] = x + x_rot;
        vertices[2*i + 1] = y + y_rot;
    } 
}

void quaternion_to_angle(double qz, double qw, double *sin_t, double *cos_t) {
    *cos_t = 2*qw*qw - 1;
    *sin_t = 2*qw*qz;
}

void rotate_point(double x, double y, double cos_t, double sin_t, double *x_rot, double *y_rot) {
    *x_rot = x*cos_t - y*sin_t;
    *y_rot = x*sin_t + y*cos_t;
}

double rand_double(double high) {
    return ((double)rand() * high) / (double)RAND_MAX;
}

int rand_int(int high) {
    return rand() % high;
}

double distance(const Particle* p1, const Particle* p2, const Grid* grid) {
    double dx = fmin(fabs(p1->x - p2->x), fabs(grid->Lx - fabs(p1->x - p2->x)));
    double dy = fmin(fabs(p1->y - p2->y), fabs(grid->Ly - fabs(p1->y - p2->y)));
    return sqrt(pow(dx, 2) + pow(dy, 2));
}

int write_pdb(const char* filename, const Grid* grid) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) return 1;
    for (size_t idx = 0; idx < grid->N; idx++)
        fprintf(file, "ATOM  %5ld  N   NONE    1    %8.3f%8.3f%8.3f  1.00  0.00\n",
            idx + 1,
            grid->points[idx].x,
            grid->points[idx].y,
            0.0
        );
    fclose(file);
    return 0;
}

// TODO: optimize write
int write_xyz(const char* filename, const Grid* grid) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) return 1;
    // 92 is the size of a line for each particle, 100 is for the header
    char *file_string = malloc((92*grid->N*(grid->n_patches + 1) + 100) * sizeof(char));
    size_t offset = 0;

    offset += sprintf(file_string + offset, "%ld\n", grid->N + grid->N * grid->n_patches);
    offset += sprintf(file_string + offset, "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3\n");
    // sprintf(file, "%ld\n", grid->N + grid->N * grid->n_patches);
    // fprintf(file, "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3\n");
    for (size_t idx = 0; idx < grid->N; idx++) {
        // fprintf(file, "%c %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
        offset += sprintf(file_string + offset, "%c %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
            grid->type,
            grid->points[idx].x,
            grid->points[idx].y,
            0.0,
            0.0,
            0.0,
            grid->points[idx].qz,
            grid->points[idx].qw,
            grid->size/2,
            grid->size/2,
            0.2
        );

        // Write patches
        for (size_t j = 0; j < grid->n_patches; j++) {
            double cos_t = 2*grid->points[idx].qw*grid->points[idx].qw - 1;
            double sin_t = 2*grid->points[idx].qw*grid->points[idx].qz;
            double x_rot, y_rot;
            rotate_point(grid->patches[j].x, grid->patches[j].y, cos_t, sin_t, &x_rot, &y_rot);
            // fprintf(file, "P %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
            offset += sprintf(file_string + offset, "P %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                grid->points[idx].x + x_rot,
                grid->points[idx].y + y_rot,
                0.0,
                0.0,
                0.0,
                grid->points[idx].qz,
                grid->points[idx].qw,
                grid->patch_size/2,
                grid->patch_size/2,
                0.2
            );
        }
    }
    fprintf(file, "%s", file_string);
    fclose(file);
    return 0;
}
