#include "simulation.h"
#include "cell.h"
#include "utils.h"
#include "movement.h"

#include <stdlib.h>
#include <stdio.h>

int simulate_random(Config *config) {
    CellLinkedGrid cll;
    random_gen(config, &cll);
    char filename[256];
    sprintf(filename, "%s/000000.xyz", config->output_folder);
    int write_code = write_xyz(filename, config, &cll);
    if (write_code) return 1;
    for (int i = 1; i <= config->num_steps; i++) {
        for (size_t idx = 0; idx < config->N; idx++) {
            size_t move_idx = rand_int(config->N);
            random_move(&cll.particles[move_idx], config, &cll);
        }

        if (i % config->save_interval == 0) {
            sprintf(filename, "%s/%06d.xyz", config->output_folder, i);
            write_code = write_xyz(filename, config, &cll);
            if (write_code) return 1;
            printf("Iteration %d finished\n", i);
        }
    }
    cll_free(&cll);
    return 0;
}
