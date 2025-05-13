#include "particles.h"
#include "utils.h"

#include <math.h>

__host__ __device__ double patch_energy(const Particle *p1, const Particle *p2, const Config* config) {
    double cos_t1 = 2*p1->qw*p1->qw - 1;
    double sin_t1 = 2*p1->qw*p1->qz;
    double cos_t2 = 2*p2->qw*p2->qw - 1;
    double sin_t2 = 2*p2->qw*p2->qz;

    for (int p1_idx = 0; p1_idx < config->num_patches; p1_idx++) {
        double x_1, y_1;
        rotate_point(config->patch_coordinates[p1_idx].x, config->patch_coordinates[p1_idx].y, cos_t1, sin_t1, &x_1, &y_1);
        x_1 += p1->x;
        y_1 += p1->y;
        Patch patch1 = {x_1, y_1};
        for (int p2_idx = 0; p2_idx < config->num_patches; p2_idx++) {
            double x_2, y_2;
            rotate_point(config->patch_coordinates[p2_idx].x, config->patch_coordinates[p2_idx].y, cos_t2, sin_t2, &x_2, &y_2);
            x_2 += p2->x;
            y_2 += p2->y;
            Patch patch2 = {x_2, y_2};
            double dx = fmin(fabs(patch1.x - patch2.x), fabs(config->Lx - fabs(patch1.x - patch2.x)));
            double dy = fmin(fabs(patch1.y - patch2.y), fabs(config->Ly - fabs(patch1.y - patch2.y)));
            if (dx * dx + dy * dy < config->patch_size * config->patch_size) {
                return config->energy_delta;
            }
        }
    }
    return 0;
}

__host__ __device__ int square_overlap(const Particle *p1, const Particle *p2, const Config *config) {
    double dx = fmin(fabs(p1->x - p2->x), fabs(config->Lx - fabs(p1->x - p2->x)));
    double dy = fmin(fabs(p1->y - p2->y), fabs(config->Ly - fabs(p1->y - p2->y)));
    if (dx * dx + dy * dy > 2.0f * config->size * config->size) return 0;
    if (dx * dx + dy * dy < config->size * config->size) return 1;

    // rotation angles from quaternions
    double cos_t1, sin_t1, cos_t2, sin_t2;
    quaternion_to_angle(p1->qz, p1->qw, &sin_t1, &cos_t1);
    quaternion_to_angle(p2->qz, p2->qw, &sin_t2, &cos_t2);

    // get vertices of squares
    double x_rot, y_rot;
    Patch patch1, patch2;
    for (int i = 0; i < 4; i++) {
        rotate_point((i % 2) ? -config->size/2 : +config->size/2, (i < 2) ? config->size/2 : -config->size/2, cos_t1, sin_t1, &x_rot, &y_rot);
        patch1 = {p1->x + x_rot, p1->y + y_rot};
        patch1.x += fabs(patch1.x - p2->x) > config->Lx/2 ? (patch1.x > p2->x) ? -config->Lx : config->Lx : 0;
        patch1.y += fabs(patch1.y - p2->y) > config->Ly/2 ? (patch1.y > p2->y) ? -config->Ly : config->Ly : 0;
        rotate_point(patch1.x - p2->x, patch1.y - p2->y, cos_t2, -sin_t2, &x_rot, &y_rot);
        if (fabs(x_rot) < config->size/2 && fabs(y_rot) < config->size/2) return 1;

        rotate_point((i % 2) ? -config->size/2 : +config->size/2, (i < 2) ? config->size/2 : -config->size/2, cos_t2, sin_t2, &x_rot, &y_rot);
        patch2 = {p2->x + x_rot, p2->y + y_rot};
        patch2.x += fabs(patch2.x - p1->x) > config->Lx/2 ? (patch2.x > p1->x) ? -config->Lx : config->Lx : 0;
        patch2.y += fabs(patch2.y - p1->y) > config->Ly/2 ? (patch2.y > p1->y) ? -config->Ly : config->Ly : 0;
        rotate_point(patch2.x - p1->x, patch2.y - p1->y, cos_t1, -sin_t1, &x_rot, &y_rot);
        if (fabs(x_rot) < config->size/2 && fabs(y_rot) < config->size/2) return 1;
    }
    
    return 0;
}

__host__ __device__ int circle_overlap(const Particle *p1, const Particle *p2, const Config *config) {
    double dx = fmin(fabs(p1->x - p2->x), fabs(config->Lx - fabs(p1->x - p2->x)));
    double dy = fmin(fabs(p1->y - p2->y), fabs(config->Ly - fabs(p1->y - p2->y)));
    double sq_size = config->size * config->size;
    return (dx * dx + dy * dy) < sq_size;
}
