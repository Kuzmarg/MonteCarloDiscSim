#ifndef MOVEMENT_H
#define MOVEMENT_H

#include "types.h"
#include <curand_kernel.h>

__host__ void random_move(Particle *p, const Config *config, CellLinkedGrid *cll);

__global__ void random_move_kernel(const Config *config, CellLinkedGrid *cll, curandState *states,
    int* shared_states, Particle* shared_particles, float* shared_energy, int stage);

#endif // MOVEMENT_H
