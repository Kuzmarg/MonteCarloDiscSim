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
    unsigned int thread_id = threadIdx.x;
    unsigned int cells_per_thread_x = cll->n_x / config->Nx_cuda, remainder_x = cll->n_x % config->Nx_cuda;
    unsigned int cells_per_thread_y = cll->n_y / config->Ny_cuda, remainder_y = cll->n_y % config->Ny_cuda;

    unsigned int x_idx = threadIdx.x / config->Ny_cuda;
    unsigned int y_idx = threadIdx.x % config->Ny_cuda;
    unsigned int x_start = x_idx * cells_per_thread_x + min(x_idx, remainder_x);
    unsigned int x_end = (x_idx + 1) * cells_per_thread_x + min(x_idx + 1, remainder_x);
    unsigned int y_start = y_idx * cells_per_thread_y + min(y_idx, remainder_y);
    unsigned int y_end = (y_idx + 1) * cells_per_thread_y + min(y_idx + 1, remainder_y);

    unsigned int X0, X1, Y0, Y1; // start and end cell coordinates for this step
    X0 = x_start + (x_end - x_start) / 2 * (stage % 2);
    Y0 = y_start + (y_end - y_start) / 2 * (stage / 2);
    X1 = (x_start + x_end) / 2 + (x_end - x_start) / 2 * (stage % 2);
    Y1 = (y_start + y_end) / 2 + (y_end - y_start) / 2 * (stage / 2);

    unsigned int rand_cell = rand_int_cuda((X1 - X0) * (Y1 - Y0), &states[thread_id]);
    unsigned int cell_x = X0 + rand_cell / (Y1 - Y0);
    unsigned int cell_y = Y0 + rand_cell % (Y1 - Y0);
    unsigned int cell_idx = cell_x * cll->n_y + cell_y;
    unsigned int cell_size = cll->head[cell_idx];
    if (cell_size == 0) return;
    unsigned int rand_particle = rand_int_cuda(cell_size, &states[thread_id]);
    Particle p = cll->cells[cell_idx * cll->max_particles + rand_particle];
    Particle moved_particle = p;

    double rand_sample = rand_double_cuda(1, &states[thread_id]);
    if (rand_sample < 0.5) {
        double angle = rand_double_cuda(2 * M_PI, &states[thread_id]);
        double distance = rand_double_cuda(config->max_translation, &states[thread_id]);
        moved_particle.x += distance * cos(angle);
        moved_particle.y += distance * sin(angle);
        moved_particle.x = fmod(moved_particle.x + config->Lx, config->Lx);
        moved_particle.y = fmod(moved_particle.y + config->Ly, config->Ly);
    } else {
        double angle = rand_double_cuda(config->max_rotation, &states[thread_id]);
        moved_particle.qw = cos(angle/2);
        moved_particle.qz = sin(angle/2);
    }
    cll_remove_point(&p, cll);
    if (cll_check_overlap(&moved_particle, cll, config)) {
        cll_add_point(&p, cll);
        return;
    }
    
    double delta_energy = cll_patch_energy(&moved_particle, cll, config) - cll_patch_energy(&p, cll, config);
    rand_sample = rand_double_cuda(1, &states[thread_id]);
    if (delta_energy > 0 && rand_sample > exp(-delta_energy)) {
        cll_add_point(&p, cll);
        return;
    }
    p = moved_particle;
    cll_add_point(&p, cll);
}
