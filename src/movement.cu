#include "movement.h"
#include "cell.h"
#include "utils.h"

#include <stdio.h>
#include <math.h>
#include <curand_kernel.h>

__host__ void random_move(Particle *p, const Config *config, CellLinkedGrid *cll) {
    double rand_sample = rand_double(1);
    Particle moved_particle = *p;
    if (rand_sample < 0.5) {
        double angle = rand_double(2 * M_PI);
        double distance = rand_double(config->max_translation);
        moved_particle.x += distance * cos(angle);
        moved_particle.y += distance * sin(angle);
        moved_particle.x = fmod(moved_particle.x + config->Lx, config->Lx);
        moved_particle.y = fmod(moved_particle.y + config->Ly, config->Ly);
    } else {
        double angle = rand_double(config->max_rotation);
        moved_particle.qw = cos(angle/2);
        moved_particle.qz = sin(angle/2);
    }

    cll_remove_point(p, cll);
    if (!cll_check_overlap(&moved_particle, cll, config)) {
        double delta_energy = cll_patch_energy(&moved_particle, cll, config) - cll_patch_energy(p, cll, config);
        rand_sample = rand_double(1);
        if (delta_energy > 0 && rand_sample > exp(-delta_energy)) {
            cll_add_point(p, cll);
            return;
        }
        *p = moved_particle;
    }
    cll_add_point(p, cll);
    return;
}

__global__ void random_move_kernel(const Config *config, CellLinkedGrid *cll, curandState* states, int stage) {
    unsigned rand_idx = threadIdx.x;
    unsigned int cells_per_thread_x = cll->n_x / config->Nx_cuda, remainder_x = cll->n_x % config->Nx_cuda;
    unsigned int cells_per_thread_y = cll->n_y / config->Ny_cuda, remainder_y = cll->n_y % config->Ny_cuda;

    unsigned int x_idx = threadIdx.x / config->Ny_cuda;
    unsigned int y_idx = threadIdx.x % config->Ny_cuda;
    unsigned int x_start = x_idx * cells_per_thread_x + min(x_idx, remainder_x);
    unsigned int x_end = (x_idx + 1) * cells_per_thread_x + min(x_idx + 1, remainder_x);
    unsigned int y_start = y_idx * cells_per_thread_y + min(y_idx, remainder_y);
    unsigned int y_end = (y_idx + 1) * cells_per_thread_y + min(y_idx + 1, remainder_y);

    unsigned int X0, X1, Y0, Y1; // start and end cell coordinates for this step
    switch (stage) {
    case 0:
        X0 = x_start;
        X1 = (x_start + x_end) / 2;
        Y0 = y_start;
        Y1 = (y_start + y_end) / 2;
        break;
    case 1:
        X0 = (x_start + x_end) / 2;
        X1 = x_end;
        Y0 = y_start;
        Y1 = (y_start + y_end) / 2;
        break;
    case 2:
        X0 = x_start;
        X1 = (x_start + x_end) / 2;
        Y0 = (y_start + y_end) / 2;
        Y1 = y_end;
        break;
    case 3:
        X0 = (x_start + x_end) / 2;
        X1 = x_end;
        Y0 = (y_start + y_end) / 2;
        Y1 = y_end;
        break;
    }

    unsigned int num_particles = 0;
    for (unsigned int x = X0; x < X1; x++)
        for (unsigned int y = Y0; y < Y1; y++)
            num_particles += cll->head[x * cll->n_y + y];
    if (num_particles == 0) return;
    unsigned int rand_particle = rand_int_cuda(num_particles, &states[rand_idx]);
    unsigned int cell_idx = 0;
    unsigned int particle_idx = 0;
    for (unsigned int x = X0; x < X1; x++) {
        for (unsigned int y = Y0; y < Y1; y++) {
            unsigned int cell_count = cll->head[x * cll->n_y + y];
            if (rand_particle < cell_count) {
                cell_idx = x * cll->n_y + y;
                particle_idx = rand_particle;
                break;
            }
            rand_particle -= cell_count;
        }
        if (rand_particle < cll->head[cell_idx]) break;
    }
    Particle p = cll->cells[cell_idx * cll->max_particles + particle_idx];
    Particle moved_particle = p;
    double rand_sample = rand_double_cuda(1, &states[rand_idx]);
    if (rand_sample < 0.5) {
        double angle = rand_double_cuda(2 * M_PI, &states[rand_idx]);
        double distance = rand_double_cuda(config->max_translation, &states[rand_idx]);
        moved_particle.x += distance * cos(angle);
        moved_particle.y += distance * sin(angle);
        moved_particle.x = fmod(moved_particle.x + config->Lx, config->Lx);
        moved_particle.y = fmod(moved_particle.y + config->Ly, config->Ly);
    } else {
        double angle = rand_double_cuda(config->max_rotation, &states[rand_idx]);
        moved_particle.qw = cos(angle/2);
        moved_particle.qz = sin(angle/2);
    }
    cll_remove_point(&p, cll);
    if (!cll_check_overlap(&moved_particle, cll, config)) {
        double delta_energy = cll_patch_energy(&moved_particle, cll, config) - cll_patch_energy(&p, cll, config);
        rand_sample = rand_double_cuda(1, &states[rand_idx]);
        if (delta_energy > 0 && rand_sample > exp(-delta_energy)) {
            cll_add_point(&p, cll);
            return;
        }
        p = moved_particle;
    }
    cll_add_point(&p, cll);
}
