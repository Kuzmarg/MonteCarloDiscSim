#include <time.h>
#include <stdlib.h>
#include <string.h>
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
