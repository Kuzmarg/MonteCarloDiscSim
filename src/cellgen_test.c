#include <time.h>
#include <stdlib.h>
#include "cell.h"

int main(int argc, char* argv[]) {
    srand(time(NULL));

    ParticleType pt = SQUARE;
    Particle *p = NULL;
    double size = 1;
    unsigned int N = 20000;
    unsigned int Lx = 200;
    unsigned int Ly = 200;
    Grid grid = {pt, p, size, N, Lx, Ly};

    random_gen("data/test.xyz", &grid);
    return 0;
}
