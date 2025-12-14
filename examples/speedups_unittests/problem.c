/**
 * Mock of shearing sheet
 *
 * This example is a simple test of collision detection
 * methods.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "rebound.h"
#include "boundary.h"
#include "mockCollision.h"
#include "tree.h"
#include "integrator_trace.h"
#include "assert.h"

static void reb_tree_get_nearest_neighbour_in_cell(struct reb_simulation* const r, struct reb_vec6d gb, struct reb_vec6d gbunmod, int ri, double p1_r, double second_largest_radius, struct reb_collision* collision_nearest, struct reb_treecell* c){
    const struct reb_particle* const particles = r->particles;
    if (c->pt>=0){
        // c is a leaf node
        int condition     = 1;
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

// This example is using a custom velocity dependend coefficient of restitution
double coefficient_of_restitution_bridges(const struct reb_simulation* const r, double v){
    // assumes v in units of [m/s]
    double eps = 0.32*pow(fabs(v)*100.,-0.234);
    if (eps>1) eps=1;
    if (eps<0) eps=0;
    return eps;
}

void heartbeat(struct reb_simulation* const r){
    if (reb_simulation_output_check(r, 1e-1*2.*M_PI/r->ri_sei.OMEGA)){
        reb_simulation_output_timing(r, 0);
        //reb_output_append_velocity_dispersion("veldisp.txt");
    }
    if (reb_simulation_output_check(r, 2.*M_PI/r->ri_sei.OMEGA)){
        //reb_simulation_output_ascii("position.txt");
    }
}

struct reb_simulation* mockShearingSheetSimulation(int collisionDetectionType, int numParticles){
    struct reb_simulation* r = reb_simulation_create();
    // Setup constants
    r->opening_angle2    = .5;                  // This determines the precission of the tree code gravity calculation.
    r->integrator        = REB_INTEGRATOR_SEI;
    r->boundary          = REB_BOUNDARY_SHEAR;
    r->gravity           = REB_GRAVITY_TREE;
    r->collision         = collisionDetectionType;
    r->collision_resolve = reb_collision_resolve_hardsphere;
    double OMEGA         = 0.00013143527;       // 1/s
    r->ri_sei.OMEGA      = OMEGA;
    r->G                 = 6.67428e-11;         // N / (1e-5 kg)^2 m^2
    r->softening         = 0.1;                 // m
    r->dt                = 1e-3*2.*M_PI/OMEGA;  // s
    r->heartbeat         = heartbeat;           // function pointer for heartbeat
    // This example uses two root boxes in the x and y direction.
    // Although not necessary in this case, it allows for the parallelization using MPI.
    // See Rein & Liu for a description of what a root box is in this context.
    double surfacedensity          = 400;     // kg/m^2
    double particle_density        = 400;     // kg/m^3
    double particle_radius_min     = 1;       // m
    double boxsize             = 100;         // m
    reb_simulation_configure_box(r, boxsize, 2, 2, 1);
    r->N_ghost_x = 2;
    r->N_ghost_y = 2;
    r->N_ghost_z = 0;

    // Use Bridges et al coefficient of restitution.
    r->coefficient_of_restitution = coefficient_of_restitution_bridges;
    // When two particles collide and the relative velocity is zero, the might sink into each other in the next time step.
    // By adding a small repulsive velocity to each collision, we prevent this from happening.
    r->minimum_collision_velocity = particle_radius_min*OMEGA*0.001;  // small fraction of the shear accross a particle


    // Add all ring paricles
    double total_mass = surfacedensity*r->boxsize.x*r->boxsize.y;
    double mass = 0;
    int i = 0;
    while(mass<total_mass && i < numParticles){
        //Make all the particles the same size and touch eachother (size 2)
        struct reb_particle pt;
        pt.x         = -r->boxsize.x/2. + i;
        pt.y         = -r->boxsize.y/2. + i;
        pt.z         = 0.5;                    // m
        pt.vx         = 0;
        pt.vy         = -1.5*pt.x*OMEGA;
        pt.vz         = 0;
        pt.ax         = 0;
        pt.ay         = 0;
        pt.az         = 0;
        double radius     = 2.0;
        pt.r         = radius;                        // m
        double        particle_mass = particle_density*4./3.*M_PI*radius*radius*radius;
        pt.m         = particle_mass;     // kg
        reb_simulation_add(r, pt);
        mass += particle_mass;
        i+= 1;
    }
    return r;
}

struct reb_simulation* mockBouncingBallsSimulation(){
    struct reb_simulation* r = reb_simulation_create();

    // Setup constants
    r->integrator    = REB_INTEGRATOR_LEAPFROG;
    r->gravity    = REB_GRAVITY_BASIC;
    r->collision    = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_hardsphere;
    r->usleep    = 1000;            // Slow down integration (for visualization only)
    r->dt = 1e-2;

    reb_simulation_configure_box(r, 3.0, 1, 1, 1);

    // Initial conditions
    {
        struct reb_particle p = {0};
        p.x  = 1; p.y  = 1; p.z  = 1;
        p.m  = 1;
        p.r  = 0.1;
        reb_simulation_add(r, p);
    }
    {
        struct reb_particle p = {0};
        p.x  = -1; p.y  = -1; p.z  = -1;
        p.m  = 1;
        p.r  = 0.1;
        reb_simulation_add(r, p);
    }
    return r;
}

double mock_tree_collisions(struct reb_simulation* r){
    r->collisions_N = 0;
    // Update and simplify tree.
    // Prepare particles for distribution to other nodes.
    reb_simulation_update_tree(r);

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
    return t1 - t0;
}

void test_Simulation(struct reb_simulation* r, double time){
    reb_simulation_integrate(r, time);
    assert(r->collisions_N == 110);
    // printf("\n\n%d\n", r->collisions_N);
}

void test_CollisionSearch(struct reb_simulation* r){
    mock_tree_collisions(r);
    //assert(r->collisions_N == 139);
    printf("\n\n%d\n", r->collisions_N);
}

void time_FullCollisionSearch(struct reb_simulation* r){
    double total_time = 0.0;
    for (int i = 0; i < 500000; i++){
        total_time += mock_reb_collision_search(r);
    }
    // assert(r->collisions_N == 110);
    printf("\nTotal Time: %f\n", total_time);
}


int main(int argc, char* argv[]){

    //Choose collision method (Case section in collision.c)
    //mockShearingSheetSimulation(Collision type (Keep REB_COLLISION_TREE, numParticles))
    struct reb_simulation* r = mockShearingSheetSimulation(REB_COLLISION_TREE, 2);
    //reb_simulation_start_server(r, 1234);
    mock_tree_collisions(r);
    test_CollisionSearch(r);
    //test_Simulation(r, 10);

    //time_CollisionSearch(r);
    reb_simulation_free(r);
}
