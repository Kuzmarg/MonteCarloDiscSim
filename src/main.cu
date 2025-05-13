#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <libconfig.h>
#include "simulation.h"

int grid_config_init_patches(Config *config, int n_patches) {
    cudaError_t result = cudaMallocManaged(&config->patch_coordinates, n_patches * sizeof(Patch));
    if (result != cudaSuccess) {
        fprintf(stderr, "Error allocating memory for patch coordinates: %s\n", cudaGetErrorString(result));
        return 1;
    }
    return 0;
}

int grid_config_destroy_patches(Config *grid) {
    return cudaFree(grid->patch_coordinates);
}

int main(int argc, char* argv[]) {
    srand(time(NULL));

    if (argc != 2) {
        printf("Usage: %s <config_file>\n", argv[0]);
        return 1;
    }

    // Read the configuration file
    config_t cfg;
    config_init(&cfg);
    if (!config_read_file(&cfg, argv[1])) {
        fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
                config_error_line(&cfg), config_error_text(&cfg));
        config_destroy(&cfg);
        return 1;
    }

    Config config;

    // CPU or CUDA run
    const char *run_type;
    if (config_lookup_string(&cfg, "run_type", &run_type)) {
        if (strcmp(run_type, "CPU") == 0) {
            config.cuda = false;
        } else if (strcmp(run_type, "CUDA") == 0) {
            config.cuda = true;
        } else {
            fprintf(stderr, "Invalid run type specified in config file\n");
            config_destroy(&cfg);
            return 1;
        }
    } else {
        fprintf(stderr, "No run type specified in config file\n");
        config_destroy(&cfg);
        return 1;
    }

    if (config.cuda)
    {
        if (config_lookup_int(&cfg, "Nx_cuda", &config.Nx_cuda) && config_lookup_int(&cfg, "Ny_cuda", &config.Ny_cuda)) {
            if (config.Nx_cuda <= 0 || config.Ny_cuda <= 0) {
                fprintf(stderr, "Invalid number of CUDA cells specified in config file\n");
                config_destroy(&cfg);
                return 1;
            }
        } else {
            fprintf(stderr, "No number of CUDA cells specified in config file\n");
            config_destroy(&cfg);
            return 1;
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////////

    const char *type_str;
    if (config_lookup_string(&cfg, "particle_type", &type_str)) {
        config.type = type_str[0] == 'C' ? CIRCLE : SQUARE;
    } else {
        fprintf(stderr, "No particle type specified in config file\n");
        config_destroy(&cfg);
        return 1;
    }

    if (config_lookup_float(&cfg, "particle_size", &config.size)) {
        if (config.size <= 0) {
            fprintf(stderr, "Invalid particle size specified in config file\n");
            config_destroy(&cfg);
            return 1;
        }
    } else {
        fprintf(stderr, "No particle size specified in config file\n");
        config_destroy(&cfg);
        return 1;
    }
    
    if (config_lookup_int(&cfg, "num_patches", &config.num_patches)) {
        if (config.num_patches < 0) {
            fprintf(stderr, "Invalid number of patches specified in config file\n");
            config_destroy(&cfg);
            return 1;
        }
    } else {
        fprintf(stderr, "No number of patches specified in config file\n");
        config_destroy(&cfg);
        return 1;
    }

    if (config_lookup_float(&cfg, "patch_size", &config.patch_size)) {
        if (config.patch_size <= 0) {
            fprintf(stderr, "Invalid patch size specified in config file\n");
            config_destroy(&cfg);
            return 1;
        }
    } else {
        fprintf(stderr, "No patch size specified in config file\n");
        config_destroy(&cfg);
        return 1;
    }
    config.patch_size *= config.size;

    if (config_lookup_float(&cfg, "energy_delta", &config.energy_delta)) {
        if (config.energy_delta > 0) {
            fprintf(stderr, "Invalid patch energy delta specified in config file\n");
            config_destroy(&cfg);
            return 1;
        }
    } else {
        fprintf(stderr, "No patch energy delta specified in config file\n");
        config_destroy(&cfg);
        return 1;
    }

    ////////////////////////////////////////////////////////////////////////////////

    if (config_lookup_int(&cfg, "num_particles", &config.N)) {
        if (config.N <= 0) {
            fprintf(stderr, "Invalid number of particles specified in config file\n");
            config_destroy(&cfg);
            return 1;
        }
    } else {
        fprintf(stderr, "No number of particles specified in config file\n");
        config_destroy(&cfg);
        return 1;
    }
    
    if (config_lookup_float(&cfg, "Lx", &config.Lx) && config_lookup_float(&cfg, "Ly", &config.Ly)) {
        if (config.Lx <= 0 || config.Ly <= 0) {
            fprintf(stderr, "Invalid box dimensions specified in config file\n");
            config_destroy(&cfg);
            return 1;
        }
    } else {
        fprintf(stderr, "No box dimensions specified in config file\n");
        config_destroy(&cfg);
        return 1;
    }

    if (config_lookup_int(&cfg, "num_steps", &config.num_steps)) {
        if (config.num_steps <= 0) {
            fprintf(stderr, "Invalid number of steps specified in config file\n");
            config_destroy(&cfg);
            return 1;
        }
    } else {
        fprintf(stderr, "No number of steps specified in config file\n");
        config_destroy(&cfg);
        return 1;
    }

    if (config_lookup_int(&cfg, "save_interval", &config.save_interval)) {
        if (config.save_interval <= 0) {
            fprintf(stderr, "Invalid save interval specified in config file\n");
            config_destroy(&cfg);
            return 1;
        }
    } else {
        fprintf(stderr, "No save interval specified in config file\n");
        config_destroy(&cfg);
        return 1;
    }

    if (config_lookup_float(&cfg, "max_rotation", &config.max_rotation)) {
        if (config.max_rotation < 0) {
            fprintf(stderr, "Invalid max rotation specified in config file\n");
            config_destroy(&cfg);
            return 1;
        }
    } else {
        fprintf(stderr, "No max rotation specified in config file\n");
        config_destroy(&cfg);
        return 1;
    }
    config.max_rotation *= M_PI;

    if (config_lookup_float(&cfg, "max_translation", &config.max_translation)) {
        if (config.max_translation < 0) {
            fprintf(stderr, "Invalid max translation specified in config file\n");
            config_destroy(&cfg);
            return 1;
        }
    } else {
        fprintf(stderr, "No max translation specified in config file\n");
        config_destroy(&cfg);
        return 1;
    }
    config.max_translation *= config.size;

    if (config_lookup_string(&cfg, "output_folder", &config.output_folder)) {
        size_t name_len = strlen(config.output_folder);
        if (name_len > 256) {
            fprintf(stderr, "Output folder name too long\n");
            config_destroy(&cfg);
            return 1;
        }
    } else {
        fprintf(stderr, "No output folder specified in config file\n");
        config_destroy(&cfg);
        return 1;
    }

    ////////////////////////////////////////////////////////////////////////////////

    if (grid_config_init_patches(&config, config.num_patches) != 0) {
        fprintf(stderr, "Failed to allocate memory for patch coordinates\n");
        config_destroy(&cfg);
        return 1;
    }
    config_setting_t *patches = config_lookup(&cfg, "patch_coordinates");
    if (patches != NULL) {
        int num_patches = config_setting_length(patches);
        if (num_patches != config.num_patches) {
            fprintf(stderr, "Number of patches in config file does not match specified number\n");
            grid_config_destroy_patches(&config);
            config_destroy(&cfg);
            return 1;
        }
        for (int i = 0; i < num_patches; i++) {
            config_setting_t *patch = config_setting_get_elem(patches, i);
            config_setting_lookup_float(patch, "x", &config.patch_coordinates[i].x);
            config.patch_coordinates[i].x *= config.size;
            config_setting_lookup_float(patch, "y", &config.patch_coordinates[i].y);
            config.patch_coordinates[i].y *= config.size;
        }
    } else {
        fprintf(stderr, "No patch coordinates specified in config file\n");
        grid_config_destroy_patches(&config);
        config_destroy(&cfg);
        return 1;
    }

    ////////////////////////////////////////////////////////////////////////////////


    printf("Running simulation with %s\n", config.cuda ? "CUDA" : "CPU");
    printf("Particle type: %s\n", config.type == CIRCLE ? "CIRCLE" : "SQUARE");
    printf("Particle size: %f\n", config.size);
    printf("Number of patches: %d\n", config.num_patches);
    printf("Patch size: %f\n", config.patch_size);
    printf("Energy delta: %f\n", config.energy_delta);
    printf("Number of particles: %d\n", config.N);
    printf("Box dimensions: %f x %f\n", config.Lx, config.Ly);
    printf("Number of steps: %d\n", config.num_steps);
    printf("Save interval: %d\n", config.save_interval);
    printf("Max rotation: %f\n", config.max_rotation);
    printf("Max translation: %f\n", config.max_translation);
    printf("Output folder: %s\n", config.output_folder);
    printf("Output energy log file: %s\n", config.output_energy_file);
    printf("Patch coordinates:\n");
    for (int i = 0; i < config.num_patches; i++) {
        printf("  Patch %d: (%f, %f)\n", i, config.patch_coordinates[i].x, config.patch_coordinates[i].y);
    }
    printf("\n");

    int result;
    if (!config.cuda)
        result = simulate_random(&config);
    else
        result = simulate_random_cuda(&config);
    grid_config_destroy_patches(&config);
    config_destroy(&cfg);
    return result;
}
