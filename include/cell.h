#ifndef CELL_GEN_H
#define CELL_GEN_H

#define MAX_GEN_ITERATIONS 10000

#include "types.h"
#include "particles.h"
#include "utils.h"

int grid_allocate(Grid *grid);

int grid_free(Grid *grid);

int cll_allocate(CellLinkedGrid *cll, const Grid *grid);

int cll_check_overlap(const Particle *p1, CellLinkedGrid *cll, const Grid* grid);

int cll_add_point(Particle *p, CellLinkedGrid *cll);

int cll_remove_point(Particle *p, CellLinkedGrid *cll);

int cll_free(CellLinkedGrid *cll);

int random_gen(Grid *grid, CellLinkedGrid *cll);

int square_gen(Grid *grid, CellLinkedGrid *cll);

int hexagonal_gen(Grid *grid, CellLinkedGrid *cll);

#endif // CELL_GEN_H
