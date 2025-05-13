#include "simulation.h"
#include "cell.h"
#include "utils.h"
#include "movement.h"

#include <stdlib.h>
#include <curand_kernel.h>
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

int simulate_random_cuda(Config *config) {
    CellLinkedGrid cll, cll_cuda;
    random_gen(config, &cll);
    cll_cuda = cll;
    cll_copy_cuda(&cll, &cll_cuda, config);
    
    Config *d_config;
    CellLinkedGrid *d_cll;
    cudaMalloc((void**)&d_config, sizeof(Config));
    cudaMemcpy(d_config, config, sizeof(Config), cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_cll, sizeof(CellLinkedGrid));
    cudaMemcpy(d_cll, &cll_cuda, sizeof(CellLinkedGrid), cudaMemcpyHostToDevice);

    // Generate initial configuration
    config->Nx_cuda = min(config->Nx_cuda, cll.n_x / 6);
    config->Ny_cuda = min(config->Ny_cuda, cll.n_y / 6);

    // Save initial configuration
    char filename[256];
    sprintf(filename, "%s/000000.xyz", config->output_folder);
    int write_code = write_xyz(filename, config, &cll);
    if (write_code) return 1;
    int n_moves = (int)ceil((double)(config->N) / (config->Nx_cuda * config->Ny_cuda * 4));

    // Initialize random number generator
    curandState *d_states;
    cudaMalloc((void**)&d_states, config->Nx_cuda * config->Ny_cuda * sizeof(curandState));
    rand_init_kernel<<<config->Ny_cuda, config->Nx_cuda>>>(d_states);
    cudaDeviceSynchronize();

    // Allocate memory for overlaps and energies
    int *shared_states;
    float *shared_energies;
    Particle *shared_particles;
    cudaMalloc((void**)&shared_states, config->Nx_cuda * config->Ny_cuda * 9 * sizeof(int));
    cudaMalloc((void**)&shared_particles, config->Nx_cuda * config->Ny_cuda * 2 * sizeof(Particle));
    cudaMalloc((void**)&shared_energies, config->Nx_cuda * config->Ny_cuda * 9 * sizeof(float));

    // Simulation steps
    for (int i = 1; i <= config->num_steps; i++) {
        for (int j = 0; j < n_moves; j++) {
            for (int stage = 0; stage < 4; stage++) {
                random_move_kernel<<<config->Ny_cuda * config->Nx_cuda, 9>>>(d_config, d_cll, d_states,
                    shared_states, shared_particles, shared_energies, stage);
                cudaDeviceSynchronize();
            }
        }

        if (i % config->save_interval == 0) {
            cudaMemcpy(&cll_cuda, d_cll, sizeof(CellLinkedGrid), cudaMemcpyDeviceToHost);
            cll_copy_host(&cll, &cll_cuda, config);
            sprintf(filename, "%s/%06d.xyz", config->output_folder, i);
            write_code = write_xyz(filename, config, &cll);
            if (write_code) return 1;
            printf("Iteration %d finished\n", i);
        }
    }
    cudaFree(shared_states);
    cudaFree(shared_particles);
    cudaFree(shared_energies);

    cll_free(&cll);
    cll_free_cuda(&cll_cuda);
    cudaFree(d_config);
    cudaFree(d_cll);
    cudaFree(d_states);
    return 0;
}
