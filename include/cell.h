#ifndef CELL_GEN_HPP
#define CELL_GEN_HPP

#include <stddef.h>
#include <stdbool.h>

typedef struct {
    double x, y;
} Point;

typedef struct {
    unsigned int N;
    double Lx, Ly;  
    double sigma;
} GenParams;

int write_pdb(const char* filename, const Point* points, size_t N);

int random_gen(const char* filename, const GenParams *params, bool pbc);

int square_gen(const char* filename, const GenParams *params, bool pbc);

int hexagonal_gen(const char* filename, const GenParams *params, bool pbc);

#endif // CELL_GEN_HPP
