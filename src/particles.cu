#include "particles.h"
#include "utils.h"

#include <math.h>

double patch_energy(const Particle *p1, const Particle *p2, const Config* config) {
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
            if (distance_patch(&patch1, &patch2, config) < config->patch_size) {
                return config->energy_delta;
            }
        }
    }
    return 0;
}

int square_overlap(const Particle *p1, const Particle *p2, const Config *config) {
    if (distance(p1, p2, config) > sqrt(2) * config->size) return 0;
    if (distance(p1, p2, config) < config->size) return 1;

    // rotation angles from quaternions
    double cos_t1, sin_t1, cos_t2, sin_t2;
    quaternion_to_angle(p1->qz, p1->qw, &sin_t1, &cos_t1);
    quaternion_to_angle(p2->qz, p2->qw, &sin_t2, &cos_t2);

    // get vertices of squares
    double vertices1[8], vertices2[8];
    get_square_vertices(p1->x, p1->y, config->size, cos_t1, sin_t1, vertices1);
    get_square_vertices(p2->x, p2->y, config->size, cos_t2, sin_t2, vertices2);

    
    // rotate vertices of square 2 by -t1 to square 1's frame of reference and check for overlap
    for (int i = 0; i < 8; i += 2) {
        double x_rot, y_rot;
        // check for boundary conditions
        Particle pv = {vertices2[i], vertices2[i + 1], 0, 0};
        if (fabs(pv.x - p1->x) > config->Lx/2) pv.x += (pv.x > p1->x) ? -config->Lx : config->Lx;
        if (fabs(pv.y - p1->y) > config->Ly/2) pv.y += (pv.y > p1->y) ? -config->Ly : config->Ly;
        rotate_point(pv.x - p1->x, pv.y - p1->y, cos_t1, -sin_t1, &x_rot, &y_rot);
        if (fabs(x_rot) < config->size/2 && fabs(y_rot) < config->size/2) return 1;
    }

    // rotate vertices of square 1 by -t2 to square 2's frame of reference and check for overlap
    for (int i = 0; i < 8; i += 2) {
        double x_rot, y_rot;
        // check for boundary conditions
        Particle pv = {vertices1[i], vertices1[i + 1], 0, 0};
        if (fabs(pv.x - p2->x) > config->Lx/2) pv.x += (pv.x > p2->x) ? -config->Lx : config->Lx;
        if (fabs(pv.y - p2->y) > config->Ly/2) pv.y += (pv.y > p2->y) ? -config->Ly : config->Ly;
        rotate_point(pv.x - p2->x, pv.y - p2->y, cos_t2, -sin_t2, &x_rot, &y_rot);
        if (fabs(x_rot) < config->size/2 && fabs(y_rot) < config->size/2) return 1;
    }
    
    return 0;
}

int circle_overlap(const Particle *p1, const Particle *p2, const Config *config) {
    return distance(p1, p2, config) < config->size;
}
