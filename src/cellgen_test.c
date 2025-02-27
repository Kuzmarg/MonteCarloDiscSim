#include <time.h>
#include <stdlib.h>
#include "cell.h"

int main(int argc, char* argv[]) {
    srand(time(NULL));

    if (argc != 7) {
        printf("Usage: %s <particle_type> <particle_size> <N> <Lx> <Ly> <output_file>\n", argv[0]);
        return 1;
    }

    ParticleType pt = atoi(argv[1]);
    Particle *p = NULL;
    double size = atof(argv[2]);
    long N = atol(argv[3]);
    double Lx = atof(argv[4]);
    double Ly = atof(argv[5]);
    Grid grid = {pt, p, size, N, Lx, Ly};

    random_gen(argv[6], &grid);
    return 0;
}
