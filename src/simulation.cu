#include "simulation.h"
#include "cell.h"
#include "utils.h"
#include "movement.h"

#include <stdlib.h>
#include <stdio.h>

int simulate_random(Grid *grid, const char *output_folder) {
    CellLinkedGrid cll;
    random_gen(grid, &cll);
    char filename[256];
    sprintf(filename, "%s/000000.xyz", output_folder);
    int write_code = write_xyz(filename, grid);
    if (write_code) return 1;
    for (int i = 0; i < grid->simulation_iterations; i++) {
        for (size_t idx = 0; idx < grid->count_move_cells; idx++) {
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
