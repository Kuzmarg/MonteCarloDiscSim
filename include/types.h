#ifndef TYPES_H
#define TYPES_H

typedef enum {
    CIRCLE = 'C',
    SQUARE = 'S',
} ParticleType;

typedef struct {
    double x, y; // position
    double qz, qw; // parts of rotation quaternion
    int cll_cell_idx; // index of cell in cll
    int id; // index of particle in grid
    double energy; // energy of the particle
} Particle;

typedef struct {
    double x, y; // position
} Patch;

typedef struct {
    ParticleType type; // type of particle
    Particle *points; // array of particles
    double size; // diameter for circles, side length for squares
    long N; // number of particles
    double Lx, Ly; // box dimensions
    Patch *patches; // array of patch coordinates that are set on each particle
    double patch_size; // size of patches
    long n_patches; // number of patches on each particle
} Grid;

typedef struct {
    Particle *cells; // array of particles of shape (n_x * n_y * max_particles)
    int *head; // array of counts of particles in each cell
    long n_x, n_y; // number of cells in x and y directions
    long max_particles; // maximum number of particles per cell
    double s_x, s_y; // size of cells in x and y directions
    int (*check_overlap)(const Particle*, const Particle*, const Grid*); // function to check overlap
} CellLinkedGrid;

#endif // TYPES_H
