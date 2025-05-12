#ifndef CELL_GEN_H
#define CELL_GEN_H

#define MAX_GEN_TRIES 10000

#include "types.h"

__host__ int cll_allocate(CellLinkedGrid *cll, const Config *config);

__host__ int cll_copy_cuda(CellLinkedGrid *cll, CellLinkedGrid *cll_cuda, Config *config);

__host__ int cll_copy_host(CellLinkedGrid *cll, CellLinkedGrid *cll_cuda, Config *config);

__host__ __device__ int cll_check_overlap(const Particle *p1, CellLinkedGrid *cll, const Config* config);

__host__ __device__ double cll_patch_energy(const Particle *p1, CellLinkedGrid *cll, const Config* config);

__host__ __device__ int cll_add_point(Particle *p, CellLinkedGrid *cll);

__host__ __device__ int cll_remove_point(Particle *p, CellLinkedGrid *cll);

__host__ int cll_free_cuda(CellLinkedGrid *cll_cuda);

__host__ int cll_free(CellLinkedGrid *cll);

__host__ int random_gen(Config *config, CellLinkedGrid *cll);

__host__ int square_gen(Config *config, CellLinkedGrid *cll);

__host__ int hexagonal_gen(Config *config, CellLinkedGrid *cll);

#endif // CELL_GEN_H
