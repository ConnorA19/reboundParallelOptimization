/**
 * @file        collision.c
 * @brief       Collision search routine.
 * @author      Hanno Rein <hanno@hanno-rein.de>
 *
 * @details     A collision is defined as an overlap between two particles. This
 * is only an approximation and works only if the timestep is small
 * enough. More precisely, dt << v / Rp, where v is the typical velocity
 * and Rp the radius of a particle. Furthermore, particles must be 
 * approaching each other at the time when they overlap. 
 * 
 * 
 * @section LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "particle.h"
#include "mockCollision.h"
#include "rebound.h"
#include "boundary.h"
#include "tree.h"
#include "integrator_trace.h"
#include <omp.h>
#ifdef MPI
#include "communication_mpi.h"
#endif // MPI
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b

static void reb_tree_get_nearest_neighbour_in_cell(struct reb_simulation* const r, struct reb_vec6d gb, struct reb_vec6d gbunmod, int ri, double p1_r,  double second_largest_radius, struct reb_collision* collision_nearest, struct reb_treecell* c);
static void reb_tree_check_for_overlapping_trajectories_in_cell(struct reb_simulation* const r, struct reb_vec6d gb, struct reb_vec6d gbunmod, int ri, double p1_r, double p1_r_plus_dtv, struct reb_collision* collision_nearest, struct reb_treecell* c, double maxdrift);

double speedup_reb_collision_search(struct reb_simulation* const r){
    double total_time = 0;
    r->collisions_N = 0;
    int N = r->N - r->N_var;
    int Ninner = N;

    int* map = NULL;
    switch (r->integrator){
        case REB_INTEGRATOR_MERCURIUS:
            if (r->ri_mercurius.mode==0){
                // After jump step, only collisions with star might occur.
                // All other collisions in encounter step/
                Ninner = 1;
            }else{
                N = r->ri_mercurius.encounter_N;
                Ninner = N;
                map = r->ri_mercurius.encounter_map;
            }
            break;
        case REB_INTEGRATOR_TRACE:
            switch (r->ri_trace.mode){
                case REB_TRACE_MODE_INTERACTION:
                case REB_TRACE_MODE_NONE:
                    // After jump step, only collisions with star might occur.
                    // All other collisions in encounter step/
                    Ninner = 1;
                    break;
                case REB_TRACE_MODE_KEPLER:
                    N = r->ri_trace.encounter_N;
                    Ninner = N;
                    map = r->ri_trace.encounter_map;
                    break;
                case REB_TRACE_MODE_FULL:
                    // Do the default collision search
                    break;
            }
            break;
        default:
            // map = NULL
            break;
    }
    const struct reb_particle* const particles = r->particles;
    switch (r->collision){
        case REB_COLLISION_NONE:
            break;
        case REB_COLLISION_DIRECT:
            {
                // Loop over ghost boxes, but only the inner most ring.
                int N_ghost_xcol = (r->N_ghost_x>1?1:r->N_ghost_x);
                int N_ghost_ycol = (r->N_ghost_y>1?1:r->N_ghost_y);
                int N_ghost_zcol = (r->N_ghost_z>1?1:r->N_ghost_z);
                for (int gbx=-N_ghost_xcol; gbx<=N_ghost_xcol; gbx++){
                    for (int gby=-N_ghost_ycol; gby<=N_ghost_ycol; gby++){
                        for (int gbz=-N_ghost_zcol; gbz<=N_ghost_zcol; gbz++){
                            // Loop over all particles
                            for (int i=0;i<N;i++){
#ifndef OPENMP
                                if (reb_sigint > 1) return;
#endif // OPENMP
                                int ip = i;
                                if (map){
                                    ip = map[i];
                                }
                                struct reb_particle p1 = particles[ip];
                                struct reb_vec6d gborig = reb_boundary_get_ghostbox(r, gbx,gby,gbz);
                                struct reb_vec6d gb = gborig;
                                // Precalculate shifted position 
                                gb.x += p1.x;
                                gb.y += p1.y;
                                gb.z += p1.z;
                                gb.vx += p1.vx;
                                gb.vy += p1.vy;
                                gb.vz += p1.vz;
                                // Loop over all particles again
                                for (int j=0;j<Ninner;j++){
                                    // Do not collide particle with itself.
                                    if (i==j) continue;
                                    int jp = j;
                                    if (map){
                                        jp = map[j];
                                    }
                                    struct reb_particle p2 = particles[jp];
                                    double dx = gb.x - p2.x; 
                                    double dy = gb.y - p2.y; 
                                    double dz = gb.z - p2.z; 
                                    double sr = p1.r + p2.r; 
                                    double r2 = dx*dx+dy*dy+dz*dz;
                                    // Check if particles are overlapping 
                                    if (r2>sr*sr) continue;    
                                    double dvx = gb.vx - p2.vx; 
                                    double dvy = gb.vy - p2.vy; 
                                    double dvz = gb.vz - p2.vz; 
                                    // Check if particles are approaching each other
                                    if (dvx*dx + dvy*dy + dvz*dz >0) continue; 
                                    // Add particles to collision array.
                                    if (r->N_allocated_collisions<=r->collisions_N){
                                        // Allocate memory if there is no space in array.
                                        // Init to 32 if no space has been allocated yet, otherwise double it.
                                        r->N_allocated_collisions = r->N_allocated_collisions ? r->N_allocated_collisions * 2 : 32;
                                        r->collisions = realloc(r->collisions,sizeof(struct reb_collision)*r->N_allocated_collisions);
                                    }
                                    r->collisions[r->collisions_N].p1 = ip;
                                    r->collisions[r->collisions_N].p2 = jp;
                                    r->collisions[r->collisions_N].gb = gborig;
                                    r->collisions_N++;
                                }
                            }
                        }
                    }
                }
            }
            break;
        case REB_COLLISION_LINE:
            {
                double dt_last_done = r->dt_last_done;
                // Loop over ghost boxes, but only the inner most ring.
                int N_ghost_xcol = (r->N_ghost_x>1?1:r->N_ghost_x);
                int N_ghost_ycol = (r->N_ghost_y>1?1:r->N_ghost_y);
                int N_ghost_zcol = (r->N_ghost_z>1?1:r->N_ghost_z);
                for (int gbx=-N_ghost_xcol; gbx<=N_ghost_xcol; gbx++){
                    for (int gby=-N_ghost_ycol; gby<=N_ghost_ycol; gby++){
                        for (int gbz=-N_ghost_zcol; gbz<=N_ghost_zcol; gbz++){
                            // Loop over all particles
                            for (int i=0;i<N;i++){
#ifndef OPENMP
                                if (reb_sigint > 1) return;
#endif // OPENMP
                                int ip = i;
                                if (map){
                                    ip = map[i];
                                }
                                struct reb_particle p1 = particles[ip];
                                struct reb_vec6d gborig = reb_boundary_get_ghostbox(r, gbx,gby,gbz);
                                struct reb_vec6d gb = gborig;
                                // Precalculate shifted position 
                                gb.x += p1.x;
                                gb.y += p1.y;
                                gb.z += p1.z;
                                gb.vx += p1.vx;
                                gb.vy += p1.vy;
                                gb.vz += p1.vz;
                                // Loop over all particles again
                                for (int j=i+1;j<N;j++){
                                    int jp = j;
                                    if (map){
                                        jp = map[j];
                                    }
                                    struct reb_particle p2 = particles[jp];
                                    const double dx1 = gb.x - p2.x; // distance at end
                                    const double dy1 = gb.y - p2.y;
                                    const double dz1 = gb.z - p2.z;
                                    const double r1 = (dx1*dx1 + dy1*dy1 + dz1*dz1);
                                    const double dvx1 = gb.vx - p2.vx; 
                                    const double dvy1 = gb.vy - p2.vy;
                                    const double dvz1 = gb.vz - p2.vz;
                                    const double dx2 = dx1 -dt_last_done*dvx1; // distance at beginning
                                    const double dy2 = dy1 -dt_last_done*dvy1;
                                    const double dz2 = dz1 -dt_last_done*dvz1;
                                    const double r2 = (dx2*dx2 + dy2*dy2 + dz2*dz2);
                                    const double t_closest = (dx1*dvx1 + dy1*dvy1 + dz1*dvz1)/(dvx1*dvx1 + dvy1*dvy1 + dvz1*dvz1);

                                    double rmin2_ab = MIN(r1,r2);
                                    if (t_closest/dt_last_done>=0. && t_closest/dt_last_done<=1.){
                                        const double dx3 = dx1-t_closest*dvx1; // closest approach
                                        const double dy3 = dy1-t_closest*dvy1;
                                        const double dz3 = dz1-t_closest*dvz1;
                                        const double r3 = (dx3*dx3 + dy3*dy3 + dz3*dz3);
                                        rmin2_ab = MIN(rmin2_ab, r3);
                                    }
                                    double rsum = p1.r + p2.r;
                                    if (rmin2_ab>rsum*rsum) continue;

                                    // Add particles to collision array.
                                    if (r->N_allocated_collisions<=r->collisions_N){
                                        // Allocate memory if there is no space in array.
                                        // Init to 32 if no space has been allocated yet, otherwise double it.
                                        r->N_allocated_collisions = r->N_allocated_collisions ? r->N_allocated_collisions * 2 : 32;
                                        r->collisions = realloc(r->collisions,sizeof(struct reb_collision)*r->N_allocated_collisions);
                                    }
                                    r->collisions[r->collisions_N].p1 = ip;
                                    r->collisions[r->collisions_N].p2 = jp;
                                    r->collisions[r->collisions_N].gb = gborig;
                                    r->collisions_N++;
                                }
                            }
                        }
                    }
                }
            }
            break;
        case REB_COLLISION_TREE:
            {
                // Update and simplify tree. 
                // Prepare particles for distribution to other nodes. 
                reb_simulation_update_tree(r);          

#ifdef MPI
                // Distribute particles and add newly received particles to tree.
                reb_communication_mpi_distribute_particles(r);

                // Prepare essential tree (and particles close to the boundary needed for collisions) for distribution to other nodes.
                reb_tree_prepare_essential_tree_for_collisions(r);

                // Transfer essential tree and particles needed for collisions.
                reb_communication_mpi_distribute_essential_tree_for_collisions(r);
#endif // MPI

                // Loop over ghost boxes, but only the inner most ring.
                int N_ghost_xcol = (r->N_ghost_x>1?1:r->N_ghost_x);
                int N_ghost_ycol = (r->N_ghost_y>1?1:r->N_ghost_y);
                int N_ghost_zcol = (r->N_ghost_z>1?1:r->N_ghost_z);
                const struct reb_particle* const particles = r->particles;
                const int N = r->N - r->N_var;
                // Find second largest radius
                int l1 = -1;
                int l2 = -1;
                reb_simulation_two_largest_particles(r, &l1, &l2);
                double second_largest_radius = 0;
                if (l2 != -1){
                    second_largest_radius = r->particles[l2].r;
                }

                //---------------------------------------------TIME HERE---------------------------
                double t0 = omp_get_wtime();
                // Loop over all particles
#pragma omp parallel for schedule(guided)
                for (int i=0;i<N;i++){
#ifndef OPENMP
                    if (reb_sigint > 1) return;
#endif // OPENMP
                    struct reb_particle p1 = particles[i];
                    struct reb_collision collision_nearest;
                    collision_nearest.p1 = i;
                    collision_nearest.p2 = -1;
                    double p1_r = p1.r;
                    // Loop over ghost boxes.
                    for (int gbx=-N_ghost_xcol; gbx<=N_ghost_xcol; gbx++){
                        for (int gby=-N_ghost_ycol; gby<=N_ghost_ycol; gby++){
                            for (int gbz=-N_ghost_zcol; gbz<=N_ghost_zcol; gbz++){
                                // Calculated shifted position (for speedup). 
                                struct reb_vec6d gb = reb_boundary_get_ghostbox(r, gbx,gby,gbz);
                                struct reb_vec6d gbunmod = gb;
                                gb.x += p1.x; 
                                gb.y += p1.y; 
                                gb.z += p1.z; 
                                gb.vx += p1.vx; 
                                gb.vy += p1.vy; 
                                gb.vz += p1.vz; 
                                // Loop over all root boxes.
                                for (int ri=0;ri<r->N_root;ri++){
                                    struct reb_treecell* rootcell = r->tree_root[ri];
                                    if (rootcell!=NULL){
                                        reb_tree_get_nearest_neighbour_in_cell(r, gb, gbunmod, ri, p1_r, second_largest_radius, &collision_nearest, rootcell);
                                    }
                                }
                            }
                        }
                    }
                    // Continue if no collision was found
                    if (collision_nearest.p2==-1) continue;
                }
                double t1 = omp_get_wtime();
                total_time = t1 - t0;
            }
            break;
        case REB_COLLISION_LINETREE:
            {
                // Calculate max drift (can also be stored in tree for further speedup)
                double vmax2 = 0.;
                for (int i=0;i<N;i++){
                    struct reb_particle p1 = particles[i];
                    vmax2 = MAX(vmax2, p1.vx*p1.vx + p1.vy*p1.vy + p1.vz*p1.vz);
                }
                double maxdrift = r->dt_last_done*sqrt(vmax2);
                // Update and simplify tree. 
                // Prepare particles for distribution to other nodes. 
                reb_simulation_update_tree(r);          

                // Loop over ghost boxes, but only the inner most ring.
                int N_ghost_xcol = (r->N_ghost_x>1?1:r->N_ghost_x);
                int N_ghost_ycol = (r->N_ghost_y>1?1:r->N_ghost_y);
                int N_ghost_zcol = (r->N_ghost_z>1?1:r->N_ghost_z);
                const struct reb_particle* const particles = r->particles;
                const int N = r->N - r->N_var;

                //-------------------------------------------TIME HERE-----------------------------
                // Loop over all particles
#pragma omp parallel for schedule(guided)
                for (int i=0;i<N;i++){
#ifndef OPENMP
                    if (reb_sigint > 1) return;
#endif // OPENMP
                    struct reb_particle p1 = particles[i];
                    struct reb_collision collision_nearest;
                    collision_nearest.p1 = i;
                    collision_nearest.p2 = -1;
                    double p1_r = p1.r;
                    // Add drift during last timestep
                    double p1_r_plus_dtv = p1_r + r->dt_last_done*sqrt(p1.vx*p1.vx + p1.vy*p1.vy + p1.vz*p1.vz);
                    // Loop over ghost boxes.
                    for (int gbx=-N_ghost_xcol; gbx<=N_ghost_xcol; gbx++){
                        for (int gby=-N_ghost_ycol; gby<=N_ghost_ycol; gby++){
                            for (int gbz=-N_ghost_zcol; gbz<=N_ghost_zcol; gbz++){
                                // Calculated shifted position (for speedup). 
                                struct reb_vec6d gb = reb_boundary_get_ghostbox(r, gbx,gby,gbz);
                                struct reb_vec6d gbunmod = gb;
                                gb.x += p1.x; 
                                gb.y += p1.y; 
                                gb.z += p1.z; 
                                gb.vx += p1.vx; 
                                gb.vy += p1.vy; 
                                gb.vz += p1.vz; 
                                // Loop over all root boxes.
                                for (int ri=0;ri<r->N_root;ri++){
                                    struct reb_treecell* rootcell = r->tree_root[ri];
                                    if (rootcell!=NULL){
                                        reb_tree_check_for_overlapping_trajectories_in_cell(r, gb, gbunmod,ri,p1_r,p1_r_plus_dtv,&collision_nearest,rootcell,maxdrift);
                                    }
                                }
                            }
                        }
                    }
                    // Continue if no collision was found
                    if (collision_nearest.p2==-1) continue;
                }
            }
            break;
        default:
            reb_exit("Collision routine not implemented.");
    }

    // randomize
    for (int i=0;i<r->collisions_N;i++){
        int new = rand_r(&(r->rand_seed))%r->collisions_N;
        struct reb_collision c1 = r->collisions[i];
        r->collisions[i] = r->collisions[new];
        r->collisions[new] = c1;
    }
    // Loop over all collisions previously found in reb_collision_search().

    enum REB_COLLISION_RESOLVE_OUTCOME (*resolve) (struct reb_simulation* const r, struct reb_collision c) = r->collision_resolve;
    if (resolve==NULL){
        // Default is to throw an exception
        resolve = reb_collision_resolve_halt;
    }
    unsigned int collision_resolve_keep_sorted = r->collision_resolve_keep_sorted;
    if (r->integrator == REB_INTEGRATOR_MERCURIUS || r->integrator == REB_INTEGRATOR_TRACE){
        collision_resolve_keep_sorted = 1; // Force keep_sorted for hybrid integrator
    }

    for (int i=0;i<r->collisions_N;i++){

        struct reb_collision c = r->collisions[i];
        if (c.p1 != -1 && c.p2 != -1){
            // Resolve collision
            enum REB_COLLISION_RESOLVE_OUTCOME outcome = resolve(r, c);

            // Remove particles
            if (outcome & REB_COLLISION_RESOLVE_OUTCOME_REMOVE_P1){
                // Remove p1
                int removedp1 = reb_simulation_remove_particle(r,c.p1,collision_resolve_keep_sorted);
                if (removedp1){
                    if (r->tree_root){ // In a tree, particles get removed later. 
                        for (int j=i+1;j<r->collisions_N;j++){ // Update other collisions
                            struct reb_collision* cp = &(r->collisions[j]);
                            // Skip collisions which involved the removed particle
                            if (cp->p1==c.p1 || cp->p2==c.p1){
                                cp->p1 = -1;
                                cp->p2 = -1;
                            }
                        }
                    }else{ // Not in a tree, particles get removed immediately 
                           // Update p2 of current collision
                        if (collision_resolve_keep_sorted){
                            if (c.p2 > c.p1){
                                c.p2--;
                            }
                        }else{
                            if (c.p2 == (int)(r->N-r->N_var)){
                                c.p2 = c.p1;
                            }
                        }
                        for (int j=i+1;j<r->collisions_N;j++){ // Update other collisions
                            struct reb_collision* cp = &(r->collisions[j]);
                            // Skip collisions which involve the removed particle
                            if (cp->p1==c.p1 || cp->p2==c.p1){
                                cp->p1 = -1;
                                cp->p2 = -1;
                            }
                            // Adjust collisions
                            if (collision_resolve_keep_sorted){
                                if (cp->p1 > c.p1){
                                    cp->p1--;
                                }
                                if (cp->p2 > c.p1){
                                    cp->p2--;
                                }
                            }else{
                                if (cp->p1 == (int)(r->N-r->N_var)){
                                    cp->p1 = c.p1;
                                }
                                if (cp->p2 == (int)(r->N-r->N_var)){
                                    cp->p2 = c.p1;
                                }
                            }
                        }
                    }
                }
            }
            if (outcome & REB_COLLISION_RESOLVE_OUTCOME_REMOVE_P2){
                // Remove p1
                int removedp2 = reb_simulation_remove_particle(r,c.p2,collision_resolve_keep_sorted);
                if (removedp2){ // Update other collisions
                    if (r->tree_root){ // In a tree, particles get removed later. 
                        for (int j=i+1;j<r->collisions_N;j++){ // Update other collisions
                            struct reb_collision* cp = &(r->collisions[j]);
                            // Skip collisions which involved the removed particle
                            if (cp->p1==c.p2 || cp->p2==c.p2){
                                cp->p1 = -1;
                                cp->p2 = -1;
                            }
                        }
                    }else{ // Not in a tree, particles get removed immediately 
                        for (int j=i+1;j<r->collisions_N;j++){
                            struct reb_collision* cp = &(r->collisions[j]);
                            // Skip collisions which involve the removed particle
                            if (cp->p1==c.p2 || cp->p2==c.p2){
                                cp->p1 = -1;
                                cp->p2 = -1;
                            }
                            // Adjust collisions
                            if (collision_resolve_keep_sorted){
                                if (cp->p1 > c.p2){
                                    cp->p1--;
                                }
                                if (cp->p2 > c.p2){
                                    cp->p2--;
                                }
                            }else{
                                if (cp->p1 == (int)(r->N-r->N_var)){
                                    cp->p1 = c.p2;
                                }
                                if (cp->p2 == (int)(r->N-r->N_var)){
                                    cp->p2 = c.p2;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return total_time;
}


/**
 * @brief Find the nearest neighbour in a cell or its daughters.
 * @details The function only returns a positive result if the particles
 * are overlapping. Thus, the name nearest neighbour is not
 * exactly true.
 * @param r REBOUND simulation to work on.
 * @param gb (Shifted) position and velocity of the particle.
 * @param ri Index of the root box currently being searched in.
 * @param p1_r Radius of the particle (this is not in gb).
 * @param second_largest_radius The radius of the second largest particles.
 * @param collision_nearest Pointer to the nearest collision found so far.
 * @param c Pointer to the cell currently being searched in.
 * @param gbunmod Ghostbox unmodified
 */
static void reb_tree_get_nearest_neighbour_in_cell(struct reb_simulation* const r, struct reb_vec6d gb, struct reb_vec6d gbunmod, int ri, double p1_r, double second_largest_radius, struct reb_collision* collision_nearest, struct reb_treecell* c){
    const struct reb_particle* const particles = r->particles;
    if (c->pt>=0){     
        // c is a leaf node
        int condition     = 1;
#ifdef MPI
        int isloc    = 1 ;
        isloc = reb_communication_mpi_rootbox_is_local(r, ri);
        if (isloc==1){
#endif // MPI
            /**
             * If this is a local cell, make sure particle is not colliding with itself.
             * If this is a remote cell, the particle number might be the same, even for 
             * different particles. 
             * TODO: This can probably be written in a cleaner way.
             */
            condition = (c->pt != collision_nearest->p1);
#ifdef MPI
        }
#endif // MPI
        if (condition){
            struct reb_particle p2;
#ifdef MPI
            if (isloc==1){
#endif // MPI
                p2 = particles[c->pt];
#ifdef MPI
            }else{
                int N_root_per_node = r->N_root/r->mpi_num;
                int proc_id = ri/N_root_per_node;
                p2 = r->particles_recv[proc_id][c->pt];
            }
#endif // MPI

            double dx = gb.x - p2.x;
            double dy = gb.y - p2.y;
            double dz = gb.z - p2.z;
            double r2 = dx*dx+dy*dy+dz*dz;
            // A closer neighbour has already been found 
            double rp = p1_r+p2.r;
            // reb_particles are not overlapping 
            if (r2 > rp*rp) return;
            double dvx = gb.vx - p2.vx;
            double dvy = gb.vy - p2.vy;
            double dvz = gb.vz - p2.vz;
            // reb_particles are not approaching each other
            if (dvx*dx + dvy*dy + dvz*dz >0) return;
            // Found a new nearest neighbour. Save it for later.
            collision_nearest->ri = ri;
            collision_nearest->p2 = c->pt;
            collision_nearest->gb = gbunmod;
            // Save collision in collisions array.
#pragma omp critical
            {
                if (r->N_allocated_collisions<=r->collisions_N){
                    // Init to 32 if no space has been allocated yet, otherwise double it.
                    r->N_allocated_collisions = r->N_allocated_collisions ? r->N_allocated_collisions * 2 : 32;
                    r->collisions = realloc(r->collisions,sizeof(struct reb_collision)*r->N_allocated_collisions);
                }
                r->collisions[r->collisions_N] = *collision_nearest;
                r->collisions_N++;
            }
        }
    }else{        
        // c is not a leaf node
        double dx = gb.x - c->x;
        double dy = gb.y - c->y;
        double dz = gb.z - c->z;
        double r2 = dx*dx + dy*dy + dz*dz;
        double rp  = p1_r + second_largest_radius + 0.86602540378443*c->w;
        // Check if we need to decent into daughter cells
        if (r2 < rp*rp ){
            for (int o=0;o<8;o++){
                struct reb_treecell* d = c->oct[o];
                if (d!=NULL){
                    reb_tree_get_nearest_neighbour_in_cell(r, gb,gbunmod,ri,p1_r,second_largest_radius,collision_nearest,d);
                }
            }
        }
    }
}


static void reb_tree_check_for_overlapping_trajectories_in_cell(struct reb_simulation* const r, struct reb_vec6d gb, struct reb_vec6d gbunmod, int ri, double p1_r, double p1_r_plus_dtv, struct reb_collision* collision_nearest, struct reb_treecell* c, double maxdrift){
    const struct reb_particle* const particles = r->particles;
    if (c->pt>=0){     
        // c is a leaf node
        if (c->pt != collision_nearest->p1){
            struct reb_particle p2 = particles[c->pt];
            double dt_done_last = r->dt_last_done;
            const double dx1 = gb.x - p2.x; // distance at beginning
            const double dy1 = gb.y - p2.y;
            const double dz1 = gb.z - p2.z;
            const double r1 = (dx1*dx1 + dy1*dy1 + dz1*dz1);
            const double dvx1 = gb.vx - p2.vx; 
            const double dvy1 = gb.vy - p2.vy;
            const double dvz1 = gb.vz - p2.vz;
            const double dx2 = dx1 -dt_done_last*dvx1; // distance at end
            const double dy2 = dy1 -dt_done_last*dvy1;
            const double dz2 = dz1 -dt_done_last*dvz1;
            const double r2 = (dx2*dx2 + dy2*dy2 + dz2*dz2);
            const double t_closest = (dx1*dvx1 + dy1*dvy1 + dz1*dvz1)/(dvx1*dvx1 + dvy1*dvy1 + dvz1*dvz1);

            double rmin2_ab = MIN(r1,r2);
            if (t_closest/dt_done_last>=0. && t_closest/dt_done_last<=1.){
                const double dx3 = dx1-t_closest*dvx1; // closest approach
                const double dy3 = dy1-t_closest*dvy1;
                const double dz3 = dz1-t_closest*dvz1;
                const double r3 = (dx3*dx3 + dy3*dy3 + dz3*dz3);
                rmin2_ab = MIN(rmin2_ab, r3);
            }
            double rsum = p1_r + p2.r;
            if (rmin2_ab>rsum*rsum) return;
            collision_nearest->ri = ri;
            collision_nearest->p2 = c->pt;
            collision_nearest->gb = gbunmod;
            // Save collision in collisions array.
#pragma omp critical
            {
                if (r->N_allocated_collisions<=r->collisions_N){
                    // Init to 32 if no space has been allocated yet, otherwise double it.
                    r->N_allocated_collisions = r->N_allocated_collisions ? r->N_allocated_collisions * 2 : 32;
                    r->collisions = realloc(r->collisions,sizeof(struct reb_collision)*r->N_allocated_collisions);
                }
                r->collisions[r->collisions_N] = *collision_nearest;
                r->collisions_N++;
            }
        }
    }else{        
        // c is not a leaf node
        double dx = gb.x - c->x;
        double dy = gb.y - c->y;
        double dz = gb.z - c->z;
        double r2 = dx*dx + dy*dy + dz*dz;
        double rp  = p1_r_plus_dtv + maxdrift + 0.86602540378443*c->w;
        // Check if we need to decent into daughter cells
        if (r2 < rp*rp ){
            for (int o=0;o<8;o++){
                struct reb_treecell* d = c->oct[o];
                if (d!=NULL){
                    reb_tree_check_for_overlapping_trajectories_in_cell(r, gb,gbunmod,ri,p1_r,p1_r_plus_dtv,collision_nearest,d,maxdrift);
                }
            }
        }
    }
}


