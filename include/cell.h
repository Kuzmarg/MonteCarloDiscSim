#ifndef CELL_GEN_H
#define CELL_GEN_H

#define MAX_GEN_TRIES 10000

#include "types.h"

int cll_allocate(CellLinkedGrid *cll, const Config *config);

int cll_check_overlap(const Particle *p1, CellLinkedGrid *cll, const Config* config);

double cll_patch_energy(const Particle *p1, CellLinkedGrid *cll, const Config* config);

int cll_add_point(Particle *p, CellLinkedGrid *cll);

int cll_remove_point(Particle *p, CellLinkedGrid *cll);

int cll_free(CellLinkedGrid *cll);

int random_gen(Config *config, CellLinkedGrid *cll);

int square_gen(Config *config, CellLinkedGrid *cll);

int hexagonal_gen(Config *config, CellLinkedGrid *cll);

#endif // CELL_GEN_H
