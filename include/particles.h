#ifndef PARTICLES_H
#define PARTICLES_H
    
#include "types.h"

double patch_energy(const Particle *p1, const Particle *p2, const Config* config);

int square_overlap(const Particle *p1, const Particle *p2, const Config *config);

int circle_overlap(const Particle *p1, const Particle *p2, const Config *config);

#endif // PARTICLES_H
