#include "cell.h"
#include "particles.h"
#include "utils.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

__host__ void set_random_particle(Particle *p, const Config *config) {
    p->x = rand_double(config->Lx);
    p->y = rand_double(config->Ly);
    if (config->type == CIRCLE) { p->qw = 1; return; }
    double angle = rand_double(M_PI / 2);
    p->qw = cos(angle/2);
    p->qz = sin(angle/2);
}

__host__ int cll_allocate(CellLinkedGrid *cll, const Config *config) {
    switch (config->type)
    {
        case CIRCLE:
            cll->s_x = config->size + config->patch_size;
            cll->s_y = config->size + config->patch_size;
            cll->max_particles = 4;
            break;
        case SQUARE:
            cll->s_x = config->size * sqrt(2) + config->patch_size;
            cll->s_y = config->size * sqrt(2) + config->patch_size;
            cll->max_particles = 5;
            break;
    }
    cll->n_x = ceil(config->Lx / cll->s_x);
    cll->n_y = ceil(config->Ly / cll->s_x);
    cll->particles = (Particle*)malloc(config->N * sizeof(Particle));
    cll->cells = (Particle*)malloc(cll->n_x * cll->n_y * cll->max_particles * sizeof(Particle));
    cll->head = (int*)malloc(cll->n_x * cll->n_y * sizeof(int));
    memset(cll->particles, 0, config->N * sizeof(Particle));
    memset(cll->cells, 0, cll->n_x * cll->n_y * cll->max_particles * sizeof(Particle));
    memset(cll->head, 0, cll->n_x * cll->n_y * sizeof(int));
    return cll->cells == NULL || cll->head == NULL || cll->particles == NULL;
}

__host__ int cll_copy_cuda(CellLinkedGrid *cll, CellLinkedGrid *cll_cuda, Config *config) {
    cudaMalloc((void**)&cll_cuda->cells, cll->n_x * cll->n_y * cll->max_particles * sizeof(Particle));
    cudaMalloc((void**)&cll_cuda->head, cll->n_x * cll->n_y * sizeof(int));
    cudaMalloc((void**)&cll_cuda->particles, config->N * sizeof(Particle));
    cudaMemcpy(cll_cuda->cells, cll->cells, cll->n_x * cll->n_y * cll->max_particles * sizeof(Particle), cudaMemcpyHostToDevice);
    cudaMemcpy(cll_cuda->head, cll->head, cll->n_x * cll->n_y * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(cll_cuda->particles, cll->particles, config->N * sizeof(Particle), cudaMemcpyHostToDevice);
    return 0;
}

__host__ int cll_copy_host(CellLinkedGrid *cll, CellLinkedGrid *cll_cuda, Config *config) {
    cudaMemcpy(cll->cells, cll_cuda->cells, cll->n_x * cll->n_y * cll->max_particles * sizeof(Particle), cudaMemcpyDeviceToHost);
    cudaMemcpy(cll->head, cll_cuda->head, cll->n_x * cll->n_y * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(cll->particles, cll_cuda->particles, config->N * sizeof(Particle), cudaMemcpyDeviceToHost);
    return 0;
}

__host__ __device__ int cll_check_overlap(const Particle *p1, CellLinkedGrid *cll, const Config* config) {
    int x_idx = (int) (p1->x / cll->s_x);
    int y_idx = (int) (p1->y / cll->s_x);
    for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
            size_t idx = (x_idx + i + cll->n_x) % cll->n_x;
            size_t idy = (y_idx + j + cll->n_y) % cll->n_y;
            size_t cell_idx = idx * cll->n_y + idy;
            for (int k = 0; k < cll->head[cell_idx]; k++) {
                switch (config->type) {
                case CIRCLE:
                    if (circle_overlap(p1, &cll->cells[cell_idx * cll->max_particles + k], config)) return 1;
                    break;
                case SQUARE:
                    if (square_overlap(p1, &cll->cells[cell_idx * cll->max_particles + k], config)) return 1;
                    break;
                }
            }
        }
    }
    return 0;
}

__host__ __device__ double cll_patch_energy(const Particle *p1, CellLinkedGrid *cll, const Config* config) {
    double energy = 0;
    int x_idx = (int) (p1->x / cll->s_x);
    int y_idx = (int) (p1->y / cll->s_x);
    for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
            size_t idx = (x_idx + i + cll->n_x) % cll->n_x;
            size_t idy = (y_idx + j + cll->n_y) % cll->n_y;
            size_t cell_idx = idx * cll->n_y + idy;
            for (int k = 0; k < cll->head[cell_idx]; k++) {
                energy += patch_energy(p1, &cll->cells[cell_idx * cll->max_particles + k], config);
            }
        }
    }
    return energy;
}

__host__ __device__ int cll_add_point(Particle *p, CellLinkedGrid *cll) {
    long x_idx = (long) (p->x / cll->s_x);
    long y_idx = (long) (p->y / cll->s_x);
    size_t cell_idx = x_idx * cll->n_y + y_idx;
    int i = cll->head[cell_idx];
    if (i >= cll->max_particles) return 1;
    p->cll_cell_idx = cell_idx;
    cll->cells[cell_idx * cll->max_particles + i] = *p;
    cll->head[cell_idx]++;
    cll->particles[p->id] = *p;
    return 0;
}

__host__ __device__ int cll_remove_point(Particle *p, CellLinkedGrid *cll) {
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

__host__ int cll_free_cuda(CellLinkedGrid *cll_cuda) {
    cudaFree(cll_cuda->cells);
    cudaFree(cll_cuda->head);
    cudaFree(cll_cuda->particles);
    return 0;
}

__host__ int cll_free(CellLinkedGrid *cll) {
    free(cll->cells);
    free(cll->head);
    free(cll->particles);
    return 0;
}

__host__ int random_gen(Config *config, CellLinkedGrid *cll) {
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
