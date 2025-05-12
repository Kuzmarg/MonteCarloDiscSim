#ifndef PARTICLES_H
#define PARTICLES_H
    
#include "types.h"

__host__ __device__ double patch_energy(const Particle *p1, const Particle *p2, const Config* config);

__host__ __device__ int square_overlap(const Particle *p1, const Particle *p2, const Config *config);

__host__ __device__ int circle_overlap(const Particle *p1, const Particle *p2, const Config *config);

#endif // PARTICLES_H
