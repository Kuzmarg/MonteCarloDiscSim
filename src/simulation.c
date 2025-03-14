#include "simulation.h"

#define SIMULATION_ITERATIONS 100
#define COUNT_MOVE_CELLS 10000

int simulate_random(Grid *grid, const char *output_folder) {
    CellLinkedGrid cll;
    random_gen(grid, &cll);
    char filename[256];
    sprintf(filename, "%s/000000.xyz", output_folder);
    int write_code = write_xyz(filename, grid);
    if (write_code) return 1;
    for (int i = 0; i < SIMULATION_ITERATIONS; i++) {
        for (size_t idx = 0; idx < COUNT_MOVE_CELLS; idx++) {
            size_t move_idx = rand_int(grid->N);
            random_move(&grid->points[move_idx], grid, &cll);
        }
        sprintf(filename, "%s/%06d.xyz", output_folder, i + 1);
        write_code = write_xyz(filename, grid);
        if (write_code) return 1;
        printf("Iteration %d finished\n", i + 1); 
    }
    grid_free(grid);
    cll_free(&cll);
    return 0;
}
