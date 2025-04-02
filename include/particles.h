#ifndef PARTICLES_H
#define PARTICLES_H
    
#include "types.h"

double patch_energy(const Particle *p1, const Particle *p2, const Grid* grid);

int square_overlap(const Particle *p1, const Particle *p2, const Grid *grid);

int circle_overlap(const Particle *p1, const Particle *p2, const Grid *grid);

#endif // PARTICLES_H
