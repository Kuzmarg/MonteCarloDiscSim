#include "cell.h"
#include "particles.h"
#include "utils.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void set_random_particle(Particle *p, const Config *grid) {
    p->x = rand_double(grid->Lx);
    p->y = rand_double(grid->Ly);
    if (grid->type == CIRCLE) { p->qw = 1; return; }
    double angle = rand_double(M_PI / 2);
    p->qw = cos(angle/2);
    p->qz = sin(angle/2);
}

int cll_allocate(CellLinkedGrid *cll, const Config *config) {
    switch (config->type)
    {
        case CIRCLE:
            cll-> s_x = config->size + config->patch_size;
            cll->max_particles = 4;
            cll->check_overlap = circle_overlap;
            break;
        case SQUARE:
            cll->s_x = config->size * sqrt(2) + config->patch_size;
            cll->max_particles = 5;
            cll->check_overlap = square_overlap;
            break;
    }
    cll->n_x = ceil(config->Lx / cll->s_x);
    cll->n_y = ceil(config->Ly / cll->s_x);
    cudaMallocManaged(&cll->particles, config->N * sizeof(Particle));
    cudaMallocManaged(&cll->cells, cll->n_x * cll->n_y * cll->max_particles * sizeof(Particle));
    cudaMallocManaged(&cll->head, cll->n_x * cll->n_y * sizeof(int));
    memset(cll->particles, 0, config->N * sizeof(Particle));
    memset(cll->cells, 0, cll->n_x * cll->n_y * cll->max_particles * sizeof(Particle));
    memset(cll->head, 0, cll->n_x * cll->n_y * sizeof(int));
    return cll->cells == NULL || cll->head == NULL || cll->particles == NULL;
}

int cll_check_overlap(const Particle *p1, CellLinkedGrid *cll, const Config* grid) {
    long x_idx = (long) (p1->x / cll->s_x);
    long y_idx = (long) (p1->y / cll->s_x);
    for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
            size_t idx = (x_idx + i + cll->n_x) % cll->n_x;
            size_t idy = (y_idx + j + cll->n_y) % cll->n_y;
            size_t cell_idx = idx * cll->n_y + idy;
            for (int k = 0; k < cll->head[cell_idx]; k++) {
                if (cll->check_overlap(p1, &cll->cells[cell_idx * cll->max_particles + k], grid)) return 1;
            }
        }
    }
    return 0;
}

double cll_patch_energy(const Particle *p1, CellLinkedGrid *cll, const Config* grid) {
    double energy = 0;
    long x_idx = (long) (p1->x / cll->s_x);
    long y_idx = (long) (p1->y / cll->s_x);
    for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
            size_t idx = (x_idx + i + cll->n_x) % cll->n_x;
            size_t idy = (y_idx + j + cll->n_y) % cll->n_y;
            size_t cell_idx = idx * cll->n_y + idy;
            for (int k = 0; k < cll->head[cell_idx]; k++) {
                energy += patch_energy(p1, &cll->cells[cell_idx * cll->max_particles + k], grid);
            }
        }
    }
    return energy;
}

int cll_add_point(Particle *p, CellLinkedGrid *cll) {
    long x_idx = (long) (p->x / cll->s_x);
    long y_idx = (long) (p->y / cll->s_x);
    size_t cell_idx = x_idx * cll->n_y + y_idx;
    int i = cll->head[cell_idx];
    if (i >= cll->max_particles) return 1;
    p->cll_cell_idx = cell_idx;
    cll->cells[cell_idx * cll->max_particles + i] = *p;
    cll->head[cell_idx]++;
    return 0;
}

int cll_remove_point(Particle *p, CellLinkedGrid *cll) {
    int cell_idx = p->cll_cell_idx;
    int cell_count = cll->head[cell_idx];
    for (int i = 0; i < cell_count; i++) {
        if (cll->cells[cell_idx * cll->max_particles + i].id == p->id) {
            cll->head[cell_idx]--;
            cll->cells[cell_idx * cll->max_particles + i] = cll->cells[cell_idx * cll->max_particles + (cell_count - 1)];
            p->cll_cell_idx = -1;
            return 0;
        }
    }
    return 1;
}

int cll_free(CellLinkedGrid *cll) {
    cudaFree(cll->cells);
    cudaFree(cll->head);
    cudaFree(cll->particles);
    return 0;
}

int random_gen(Config *config, CellLinkedGrid *cll) {
    if (cll_allocate(cll, config) != 0) return 1;

    size_t num_generated = config->N;
    for (size_t idx = 0; idx < config->N; idx++) {
        cll->particles[idx].id = idx;
        for (size_t i = 0; i < MAX_GEN_TRIES; i++) {
            set_random_particle(&cll->particles[idx], config);
            if (!cll_check_overlap(&cll->particles[idx], cll, config)) {
                cll_add_point(&cll->particles[idx], cll);
                break;
            }
            if (i == MAX_GEN_TRIES - 1) num_generated = idx;
        }
        if (num_generated < config->N) break;
    }
    printf("Generated %lu particles\n", num_generated);
    config->N = num_generated;
    return 0;
}

// int square_gen(Config *grid, CellLinkedGrid *cll) {
//     grid_allocate(grid);
//     long N_x = (long) (grid->Lx / grid->size) + 1;
//     long N_y = (long) (grid->Ly / grid->size) + 1;
//     grid->N = N_x * N_y;
//     for (long x_idx = 0; x_idx < N_x; x_idx++) {
//         for (long y_idx = 0; y_idx < N_y; y_idx++) {
//             grid->points[x_idx * N_y + y_idx].x = x_idx * grid->size;
//             grid->points[x_idx * N_y + y_idx].y = y_idx * grid->size;
//             grid->points[x_idx * N_y + y_idx].qw = 1;
//         }
//     }
//     return 0;
// }

// int hexagonal_gen(Config *grid, CellLinkedGrid *cll) {
//     long N_x = (long) ceil(grid->Lx / grid->size / (sqrt(3)/2));
//     long N_y = (long) (grid->Ly / grid->size) + 1;
//     grid-> N = N_x * N_y;
//     grid_allocate(grid);
//     for (long x_idx = 0; x_idx < N_x; x_idx++) {
//         for (long y_idx = 0; y_idx < N_y; y_idx++) {
//             grid->points[x_idx * N_y + y_idx].x = x_idx * grid->size * sqrt(3)/2;
//             grid->points[x_idx * N_y + y_idx].y = y_idx * grid->size;
//             grid->points[x_idx * N_y + y_idx].qw = 1;
//             if (x_idx % 2) grid->points[x_idx * N_y + y_idx].y += grid->size / 2;
//         }
//     }
//     return 0;
// }
