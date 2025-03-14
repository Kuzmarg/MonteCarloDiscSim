#ifndef TYPES_H
#define TYPES_H

typedef enum {
    CIRCLE = 'C',
    SQUARE = 'S',
} ParticleType;

struct Particle {
    double x, y; // position
    double qz, qw; // parts of rotation quaternion
    struct Particle *cll_copy; // pointer to particle in cell linked list
};
typedef struct Particle Particle;

typedef struct {
    ParticleType type; // type of particle
    Particle *points; // array of particles
    double size; // diameter for circles, side length for squares
    long N; // number of particles
    double Lx, Ly; // box dimensions
} Grid;

typedef struct {
    Particle *cells; // array of particles of shape (n_x * n_y * max_particles)
    long n_x, n_y; // number of cells in x and y directions
    long max_particles; // maximum number of particles per cell
    double s_x, s_y; // size of cells in x and y directions
    int (*check_overlap)(const Particle*, const Particle*, const Grid*); // function to check overlap
} CellLinkedGrid;

#endif // TYPES_H
