#include "cell.h"
#include <time.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    const unsigned int N = 380;
    const double Lx = 24.0, Ly = 24.0;
    const double sigma = 1.0;
    const GenParams params = {N, Lx, Ly, sigma};
    
    srand(time(NULL));
    random_gen("random_gen.pdb", &params, true);
    square_gen("square_gen.pdb", &params, true);
    hexagonal_gen("hexagonal_gen.pdb", &params, true);
    return 0;
}
