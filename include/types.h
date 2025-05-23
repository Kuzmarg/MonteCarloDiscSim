#ifndef TYPES_H
#define TYPES_H

typedef enum {
    CIRCLE,
    SQUARE,
} ParticleType;

typedef struct {
    double x, y; // position
    double qz, qw; // parts of rotation quaternion
    int cll_cell_idx; // index of cell in cll
    int id; // index of particle in grid
} Particle;

typedef struct {
    double x, y; // position
} Patch;

typedef struct {
    bool cuda; // true if using CUDA

    ParticleType type; // type of particle
    double size; // diameter for circles, side length for squares
    int num_patches; // number of patches on the particle
    double patch_size; // size of patches
    double energy_delta; // energy of the patch interaction
    Patch *patch_coordinates; // array of patch coordinates that are set on each particle

    double max_rotation; // maximum rotation of the particle per step
    double max_translation; // maximum translation of the particle per step

    int N; // number of particles
    double Lx, Ly; // box dimensions
    int num_steps; // number of iterations for the simulation
    int save_interval; // interval for saving the simulation
    int Nx_cuda, Ny_cuda; // number of cells in x and y directions for CUDA

    const char *output_folder; // folder to save the simulation
    const char *output_energy_file; // file to save the energy
} Config;

typedef struct {
    Particle *particles; // array of all particles in the simulation
    Particle *cells; // array of particles of shape (n_x * n_y * max_particles)
    int *head; // array of counts of particles in each cell
    int n_x, n_y; // number of cells in x and y directions
    int max_particles; // maximum number of particles per cell
    double s_x, s_y; // size of cells in x and y directions
} CellLinkedGrid;

#endif // TYPES_H
