#ifndef MOVEMENT_H
#define MOVEMENT_H

#include "types.h"
#include "particles.h"
#include "utils.h"
#include "cell.h"

int random_move(Particle *p, const Grid *grid, CellLinkedGrid *cll);

#endif // MOVEMENT_H
