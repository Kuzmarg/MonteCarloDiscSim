#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "simulation.h"

int main(int argc, char* argv[]) {
    srand(time(NULL));

    if (argc != 7) {
        printf("Usage: %s <particle_type> <particle_size> <N> <Lx> <Ly> <output_folder>\n", argv[0]);
        return 1;
    }

    ParticleType pt = atoi(argv[1]);
    Particle *p = NULL;
    double size = atof(argv[2]);
    long N = atol(argv[3]);
    double Lx = atof(argv[4]);
    double Ly = atof(argv[5]);
    Grid grid = {pt, p, size, N, Lx, Ly};

    // For testing of patches
    grid.n_patches = 4;
    grid.patch_size = size / 5;
    grid.delta_energy = 0.1;
    grid.patches = malloc(grid.n_patches * sizeof(Patch));
    grid.patches[0] = (Patch){size/2, 0};
    grid.patches[1] = (Patch){0, size/2};
    grid.patches[2] = (Patch){-size/2, 0};
    grid.patches[3] = (Patch){0, -size/2};

    char *output_folder = argv[6];
    size_t name_len = strlen(output_folder);
    
    if (name_len > 256) {
        printf("Output folder name too long\n");
        return 1;
    }
    if (output_folder[name_len - 1] == '/') {
        output_folder[name_len - 1] = '\0';
    }
    
    simulate_random(&grid, output_folder);
    return 0;
}
