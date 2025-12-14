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
#include "speedupCollisionAttempts.h"
#include "tree.h"
#include "integrator_trace.h"
#include "assert.h"

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


//Tests  ---------------------------------


void test_Simulation(){
    struct reb_simulation* r = mockShearingSheetSimulation(REB_COLLISION_TREE, 2);
    reb_simulation_integrate(r, 10);
    assert(r->collisions_N == 110);
    // printf("\n\n%d\n", r->collisions_N);
}

void test_CollisionSearchOriginal(){
    struct reb_simulation* r = mockShearingSheetSimulation(REB_COLLISION_TREE, 8);
    mock_reb_collision_search(r);
    // assert(r->collisions_N == 34);
    printf("\n\n%d\n", r->collisions_N);
}


void test_CollisionSearchCompare(int particleCount){
    struct reb_simulation* r = mockShearingSheetSimulation(REB_COLLISION_TREE, particleCount);
    mock_reb_collision_search(r);
    int collisions_original = r->collisions_N;
    reb_simulation_free(r);
    r = mockShearingSheetSimulation(REB_COLLISION_TREE, particleCount);
    speedup_reb_collision_search(r);
    int collisions_new = r->collisions_N;
    assert(collisions_original == collisions_new);
    reb_simulation_free(r);
}


void time_CollisionSearchOriginal(){
    double total_time = 0.0;
    struct reb_simulation* r;
    for (int i = 0; i < 50000; i++){
        r = mockShearingSheetSimulation(REB_COLLISION_TREE, i % 200);
        total_time += mock_reb_collision_search(r);
        reb_simulation_free(r);
    }
    // assert(r->collisions_N == 110);
    printf("\nTotal Time: %f\n", total_time);
}

void time_CollisionSearchSpeedup(){
    double total_time = 0.0;
    struct reb_simulation* r;
    for (int i = 0; i < 50000; i++){
        r = mockShearingSheetSimulation(REB_COLLISION_TREE, i % 200);
        total_time += speedup_reb_collision_search(r);
        reb_simulation_free(r);
    }
    // assert(r->collisions_N == 110);
    printf("\nTotal Time: %f\n", total_time);
}

int main(int argc, char* argv[]){

    //Choose collision method (Case section in collision.c)
    //mockShearingSheetSimulation(Collision type (Keep REB_COLLISION_TREE, numParticles))
    //reb_simulation_start_server();
    // mock_tree_collisions();

    //TEST
    test_CollisionSearchCompare(-1);
    test_CollisionSearchCompare(0);
    test_CollisionSearchCompare(1);
    test_CollisionSearchCompare(2);
    test_CollisionSearchCompare(10);
    test_CollisionSearchCompare(100);
    test_CollisionSearchCompare(1000);
    //WANT TO ADD END TO END TEST HERE WITH OUTPUT TO BIN

    //TIME
    time_CollisionSearchSpeedup();


    //test_Simulation(r, 10);

    // struct reb_simulation* r = mockShearingSheetSimulation(REB_COLLISION_TREE, 8);
    // time_FullCollisionSearchOriginal(r);
    //time_CollisionSearchOriginal(r);
    // reb_simulation_free(r);

}
