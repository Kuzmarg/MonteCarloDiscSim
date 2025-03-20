#include "cell.h"

void set_random_particle(Particle *p, const Grid *grid) {
    p->x = rand_double(grid->Lx);
    p->y = rand_double(grid->Ly);
    if (grid->type == CIRCLE) { p->qw = 1; return; }
    double angle = rand_double(M_PI / 2);
    p->qw = cos(angle/2);
    p->qz = sin(angle/2);
}

int grid_allocate(Grid *grid) {
    grid->points = (Particle*)calloc(grid->N, sizeof(Particle));
    return grid->points == NULL;
}

int grid_free(Grid *grid) {
    free(grid->points);
    return 0;
}

int cll_allocate(CellLinkedGrid *cll, const Grid *grid) {
    switch (grid->type)
    {
        case CIRCLE:
            cll-> s_x = grid->size;
            cll->max_particles = 4;
            cll->check_overlap = circle_overlap;
            break;
        case SQUARE:
            cll->s_x = grid->size * sqrt(2);
            cll->max_particles = 5;
            cll->check_overlap = square_overlap;
            break;
    }
    cll->n_x = ceil(grid->Lx / cll->s_x);
    cll->n_y = ceil(grid->Ly / cll->s_x);
    cll->cells = (Particle*)calloc(cll->n_x * cll->n_y * cll->max_particles, sizeof(Particle));
    cll->head = (int*)calloc(cll->n_x * cll->n_y, sizeof(int));
    return cll->cells == NULL;
}

int cll_check_overlap(const Particle *p1, CellLinkedGrid *cll, const Grid* grid) {
    long x_idx = (long) (p1->x / cll->s_x);
    long y_idx = (long) (p1->y / cll->s_x);
    for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
            size_t idx = (x_idx + i + cll->n_x) % cll->n_x;
            size_t idy = (y_idx + j + cll->n_y) % cll->n_y;
            size_t cell_idx = idx * cll->n_y + idy;
            for (int k = 0; k < cll->head[cell_idx]; k++) {
                if (cll->check_overlap(p1, &cll->cells[cell_idx * cll->max_particles + k], grid)) return 1;
            }
        }
    }
    return 0;
}

int cll_add_point(Particle *p, CellLinkedGrid *cll) {
    long x_idx = (long) (p->x / cll->s_x);
    long y_idx = (long) (p->y / cll->s_x);
    size_t cell_idx = x_idx * cll->n_y + y_idx;
    int i = cll->head[cell_idx];
    if (i >= cll->max_particles) return 1;
    p->cll_cell_idx = cell_idx;
    cll->cells[cell_idx * cll->max_particles + i] = *p;
    cll->head[cell_idx] = (i + 1);
    return 0;
}

int cll_remove_point(Particle *p, CellLinkedGrid *cll) {
    int cell_idx = p->cll_cell_idx;
    int cell_count = cll->head[cell_idx];
    for (int i = 0; i < cell_count; i++) {
        if (cll->cells[cell_idx * cll->max_particles + i].id == p->id) {
            cll->head[cell_idx]--;
            cll->cells[cell_idx * cll->max_particles + i] = cll->cells[cell_idx * cll->max_particles + (cell_count - 1)];
            p->cll_cell_idx = -1;
            return 0;
        }
    }
    return 1;
}

int cll_free(CellLinkedGrid *cll) {
    free(cll->cells);
    free(cll->head);
    return 0;
}

int random_gen(Grid *grid, CellLinkedGrid *cll) {
    cll_allocate(cll, grid);
    if (!cll->cells) return 1;
    grid_allocate(grid);
    if (!grid->points) { cll_free(cll); return 1; }

    size_t num_generated = grid->N;
    for (size_t idx = 0; idx < grid->N; idx++) {
        grid->points[idx].id = idx;
        for (size_t i = 0; i < MAX_GEN_ITERATIONS; i++) {
            set_random_particle(&grid->points[idx], grid);
            if (!cll_check_overlap(&grid->points[idx], cll, grid)) {
                cll_add_point(&grid->points[idx], cll);
                break;
            }
            if (i == MAX_GEN_ITERATIONS - 1) num_generated = idx;
        }
        if (num_generated < grid->N) break;
    }
    printf("Generated %lu particles\n", num_generated);
    grid->N = num_generated;
    return 0;
}

int square_gen(Grid *grid, CellLinkedGrid *cll) {
    grid_allocate(grid);
    long N_x = (long) (grid->Lx / grid->size) + 1;
    long N_y = (long) (grid->Ly / grid->size) + 1;
    grid->N = N_x * N_y;
    for (long x_idx = 0; x_idx < N_x; x_idx++) {
        for (long y_idx = 0; y_idx < N_y; y_idx++) {
            grid->points[x_idx * N_y + y_idx].x = x_idx * grid->size;
            grid->points[x_idx * N_y + y_idx].y = y_idx * grid->size;
            grid->points[x_idx * N_y + y_idx].qw = 1;
        }
    }
    return 0;
}

// TODO: fix slight boundary overstep
int hexagonal_gen(Grid *grid, CellLinkedGrid *cll) {
    long N_x = (long) ceil(grid->Lx / grid->size / (sqrt(3)/2));
    long N_y = (long) (grid->Ly / grid->size) + 1;
    grid-> N = N_x * N_y;
    grid_allocate(grid);
    for (long x_idx = 0; x_idx < N_x; x_idx++) {
        for (long y_idx = 0; y_idx < N_y; y_idx++) {
            grid->points[x_idx * N_y + y_idx].x = x_idx * grid->size * sqrt(3)/2;
            grid->points[x_idx * N_y + y_idx].y = y_idx * grid->size;
            grid->points[x_idx * N_y + y_idx].qw = 1;
            if (x_idx % 2) grid->points[x_idx * N_y + y_idx].y += grid->size / 2;
        }
    }
    return 0;
}
