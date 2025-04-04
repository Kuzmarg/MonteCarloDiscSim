#include "movement.h"
#include "cell.h"
#include "utils.h"

#include <math.h>

int random_move(Particle *p, const Grid *grid, CellLinkedGrid *cll) {
    double rand_sample = rand_double(1);
    Particle moved_particle = *p;
    if (rand_sample < 0.5) {
        double angle = rand_double(2 * M_PI);
        double distance = rand_double(grid->size);
        moved_particle.x += distance * cos(angle);
        moved_particle.y += distance * sin(angle);
    } else {
        double angle = rand_double(M_PI/2);
        moved_particle.qw = cos(angle/2);
        moved_particle.qz = sin(angle/2);
    }

    cll_remove_point(p, cll);
    moved_particle.x = fmod(moved_particle.x + grid->Lx, grid->Lx);
    moved_particle.y = fmod(moved_particle.y + grid->Ly, grid->Ly);
    if (!cll_check_overlap(&moved_particle, cll, grid)) {
        double delta_energy = cll_patch_energy(&moved_particle, cll, grid) - cll_patch_energy(p, cll, grid);
        rand_sample = rand_double(1);
        if (delta_energy > 0 && rand_sample > exp(-delta_energy)) {
            cll_add_point(p, cll);
            return 0;
        }
        *p = moved_particle;
    }
    cll_add_point(p, cll);
    return 0;
}
