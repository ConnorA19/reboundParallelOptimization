/**
 * @file collision.c
 * @brief Tree-based collision detection (Serial + OpenMP + MPI safe)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "rebound.h"
#include "tree.h"
#include "boundary.h"

#ifdef MPI
#include "communication_mpi.h"
#endif

/* ============================================================
 * Internal helper: thread-safe recursive traversal
 * ============================================================ */
static void reb_tree_find_collisions_in_cell_omp(
    struct reb_simulation* r,
    const struct reb_treecell* c,
    int ri,
    struct reb_collision** local_collisions,
    int* local_counts,
    int* local_capacities
){
    if (!c) return;

#ifdef _OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif

    /* ---------------- Leaf ---------------- */
    if (c->pt >= 0){
        if (c->pt == ri) return;

        struct reb_particle p1 = r->particles[ri];
        struct reb_particle p2;

#ifdef MPI
        if (reb_communication_mpi_rootbox_is_local(r, ri)){
            p2 = r->particles[c->pt];
        } else {
            int N_root_per_node = r->N_root / r->mpi_num;
            int proc_id = ri / N_root_per_node;
            if (proc_id >= r->mpi_num) proc_id = r->mpi_num - 1;
            p2 = r->particles_recv[proc_id][c->pt];
        }
#else
        p2 = r->particles[c->pt];
#endif

        const double dx = p1.x - p2.x;
        const double dy = p1.y - p2.y;
        const double dz = p1.z - p2.z;
        const double rp = p1.r + p2.r;

        if (dx*dx + dy*dy + dz*dz > rp*rp) return;

        const double dvx = p1.vx - p2.vx;
        const double dvy = p1.vy - p2.vy;
        const double dvz = p1.vz - p2.vz;

        if (dx*dvx + dy*dvy + dz*dvz >= 0.) return;

        /* Grow buffer */
        if (local_counts[tid] == local_capacities[tid]){
            local_capacities[tid] *= 2;
            local_collisions[tid] = realloc(
                local_collisions[tid],
                local_capacities[tid] * sizeof(**local_collisions)
            );
        }

        struct reb_collision col;
        col.p1 = ri;
        col.p2 = c->pt;
        col.ri = ri;

        local_collisions[tid][local_counts[tid]++] = col;
        return;
    }

    /* ---------------- Internal node ---------------- */
    for (int k = 0; k < 8; k++){
        if (c->oct[k]){
            reb_tree_find_collisions_in_cell_omp(
                r, c->oct[k], ri,
                local_collisions,
                local_counts,
                local_capacities
            );
        }
    }
}

/* ============================================================
 * Public entry point
 * ============================================================ */
void reb_tree_check_for_collisions(struct reb_simulation* r){
    struct reb_collision** local_collisions = NULL;
    int* local_counts = NULL;
    int* local_capacities = NULL;
    int nthreads = 1;

#ifdef _OPENMP
    /* Allocate thread-local buffers correctly */
    #pragma omp parallel
    {
        #pragma omp single
        {
            nthreads = omp_get_num_threads();

            local_collisions = calloc(nthreads, sizeof(*local_collisions));
            local_counts = calloc(nthreads, sizeof(*local_counts));
            local_capacities = calloc(nthreads, sizeof(*local_capacities));

            for (int t = 0; t < nthreads; t++){
                local_capacities[t] = 16;
                local_counts[t] = 0;
                local_collisions[t] = malloc(
                    local_capacities[t] * sizeof(**local_collisions)
                );
            }
        }
    }

    #pragma omp parallel for schedule(guided)
    for (int ri = 0; ri < r->N_root; ri++){
        reb_tree_find_collisions_in_cell_omp(
            r,
            r->root,
            ri,
            local_collisions,
            local_counts,
            local_capacities
        );
    }
#else
    /* Serial fallback */
    local_collisions = calloc(1, sizeof(*local_collisions));
    local_counts = calloc(1, sizeof(*local_counts));
    local_capacities = calloc(1, sizeof(*local_capacities));

    local_capacities[0] = 16;
    local_counts[0] = 0;
    local_collisions[0] = malloc(
        local_capacities[0] * sizeof(**local_collisions)
    );

    for (int ri = 0; ri < r->N_root; ri++){
        reb_tree_find_collisions_in_cell_omp(
            r,
            r->root,
            ri,
            local_collisions,
            local_counts,
            local_capacities
        );
    }
#endif

    /* ---------------- Merge results ---------------- */
    int total = 0;
    for (int t = 0; t < nthreads; t++){
        total += local_counts[t];
    }

    r->collisions = realloc(
        r->collisions,
        (r->collisions_N + total) * sizeof(*r->collisions)
    );

    int offset = r->collisions_N;
    for (int t = 0; t < nthreads; t++){
        memcpy(&r->collisions[offset],
               local_collisions[t],
               local_counts[t] * sizeof(*r->collisions));
        offset += local_counts[t];
    }
    r->collisions_N += total;

    /* ---------------- Cleanup ---------------- */
    for (int t = 0; t < nthreads; t++){
        free(local_collisions[t]);
    }
    free(local_collisions);
    free(local_counts);
    free(local_capacities);
}
