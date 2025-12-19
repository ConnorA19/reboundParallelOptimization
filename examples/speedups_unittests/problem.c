/**
 * Mock of shearing sheet
 *
 * This example is a simple test of collision detection
 * methods.
 */
#include "assert.h"
#include "boundary.h"
#include "integrator_trace.h"
#include "mockCollision.h"
#include "rebound.h"
#include "speedupCollisionAttempts.h"
#include "tree.h"
#include <fcntl.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// This example is using a custom velocity dependend coefficient of restitution
double coefficient_of_restitution_bridges(const struct reb_simulation *const r,
                                          double v) {
  // assumes v in units of [m/s]
  double eps = 0.32 * pow(fabs(v) * 100., -0.234);
  if (eps > 1)
    eps = 1;
  if (eps < 0)
    eps = 0;
  return eps;
}

void heartbeat(struct reb_simulation *const r) {
  if (reb_simulation_output_check(r, 1e-1 * 2. * M_PI / r->ri_sei.OMEGA)) {
    reb_simulation_output_timing(r, 0);
    // reb_output_append_velocity_dispersion("veldisp.txt");
  }
  if (reb_simulation_output_check(r, 2. * M_PI / r->ri_sei.OMEGA)) {
    // reb_simulation_output_ascii("position.txt");
  }
}

struct reb_simulation *mockShearingSheetSimulation(int collisionDetectionType,
                                                   int numParticles) {
  struct reb_simulation *r = reb_simulation_create();
  // Setup constants
  r->opening_angle2 = .5; // This determines the precission of the tree code
                          // gravity calculation.
  r->integrator = REB_INTEGRATOR_SEI;
  r->boundary = REB_BOUNDARY_SHEAR;
  r->gravity = REB_GRAVITY_TREE;
  r->collision = collisionDetectionType;
  r->collision_resolve = reb_collision_resolve_hardsphere;
  double OMEGA = 0.00013143527; // 1/s
  r->ri_sei.OMEGA = OMEGA;
  r->G = 6.67428e-11;               // N / (1e-5 kg)^2 m^2
  r->softening = 0.1;               // m
  r->dt = 1e-3 * 2. * M_PI / OMEGA; // s
  r->heartbeat = heartbeat;         // function pointer for heartbeat
  // This example uses two root boxes in the x and y direction.
  // Although not necessary in this case, it allows for the parallelization
  // using MPI. See Rein & Liu for a description of what a root box is in this
  // context.
  double surfacedensity = 400;    // kg/m^2
  double particle_density = 400;  // kg/m^3
  double particle_radius_min = 1; // m
  double boxsize = 100;           // m
  reb_simulation_configure_box(r, boxsize, 2, 2, 1);
  r->N_ghost_x = 2;
  r->N_ghost_y = 2;
  r->N_ghost_z = 0;

  // Use Bridges et al coefficient of restitution.
  r->coefficient_of_restitution = coefficient_of_restitution_bridges;
  // When two particles collide and the relative velocity is zero, the might
  // sink into each other in the next time step. By adding a small repulsive
  // velocity to each collision, we prevent this from happening.
  r->minimum_collision_velocity =
      particle_radius_min * OMEGA *
      0.001; // small fraction of the shear accross a particle

  // Add all ring paricles
  double total_mass = surfacedensity * r->boxsize.x * r->boxsize.y;
  double mass = 0;
  int i = 0;
  while (mass < total_mass && i < numParticles) {
    // Make all the particles the same size and touch eachother (size 2)
    struct reb_particle pt;
    pt.x = -r->boxsize.x / 2. + i;
    pt.y = -r->boxsize.y / 2. + i;
    pt.z = 0.5; // m
    pt.vx = 0;
    pt.vy = -1.5 * pt.x * OMEGA;
    pt.vz = 0;
    pt.ax = 0;
    pt.ay = 0;
    pt.az = 0;
    double radius = 2.0;
    pt.r = radius; // m
    double particle_mass =
        particle_density * 4. / 3. * M_PI * radius * radius * radius;
    pt.m = particle_mass; // kg
    reb_simulation_add(r, pt);
    mass += particle_mass;
    i += 1;
  }
  return r;
}

struct reb_simulation *mockBouncingBallsSimulation() {
  struct reb_simulation *r = reb_simulation_create();

  // Setup constants
  r->integrator = REB_INTEGRATOR_LEAPFROG;
  r->gravity = REB_GRAVITY_BASIC;
  r->collision = REB_COLLISION_DIRECT;
  r->collision_resolve = reb_collision_resolve_hardsphere;
  r->usleep = 1000; // Slow down integration (for visualization only)
  r->dt = 1e-2;

  reb_simulation_configure_box(r, 3.0, 1, 1, 1);

  // Initial conditions
  {
    struct reb_particle p = {0};
    p.x = 1;
    p.y = 1;
    p.z = 1;
    p.m = 1;
    p.r = 0.1;
    reb_simulation_add(r, p);
  }
  {
    struct reb_particle p = {0};
    p.x = -1;
    p.y = -1;
    p.z = -1;
    p.m = 1;
    p.r = 0.1;
    reb_simulation_add(r, p);
  }
  return r;
}

struct reb_simulation *
mockUnbalancedShearingSheetSimulation(int collisionDetectionType,
                                      int numParticles) {
  struct reb_simulation *r = reb_simulation_create();
  // Setup constants
  r->opening_angle2 = .5; // This determines the precission of the tree code
                          // gravity calculation.
  r->integrator = REB_INTEGRATOR_SEI;
  r->boundary = REB_BOUNDARY_SHEAR;
  r->gravity = REB_GRAVITY_TREE;
  r->collision = collisionDetectionType;
  r->collision_resolve = reb_collision_resolve_hardsphere;
  double OMEGA = 0.00013143527; // 1/s
  r->ri_sei.OMEGA = OMEGA;
  r->G = 6.67428e-11;               // N / (1e-5 kg)^2 m^2
  r->softening = 0.1;               // m
  r->dt = 1e-3 * 2. * M_PI / OMEGA; // s
  r->heartbeat = heartbeat;         // function pointer for heartbeat
  // This example uses two root boxes in the x and y direction.
  // Although not necessary in this case, it allows for the parallelization
  // using MPI. See Rein & Liu for a description of what a root box is in this
  // context.
  double surfacedensity = 400;    // kg/m^2
  double particle_density = 400;  // kg/m^3
  double particle_radius_min = 1; // m
  double boxsize = 100;           // m
  reb_simulation_configure_box(r, boxsize, 2, 2, 1);
  r->N_ghost_x = 2;
  r->N_ghost_y = 2;
  r->N_ghost_z = 0;

  // Use Bridges et al coefficient of restitution.
  r->coefficient_of_restitution = coefficient_of_restitution_bridges;
  // When two particles collide and the relative velocity is zero, the might
  // sink into each other in the next time step. By adding a small repulsive
  // velocity to each collision, we prevent this from happening.
  r->minimum_collision_velocity =
      particle_radius_min * OMEGA *
      0.001; // small fraction of the shear accross a particle

  // Add all ring paricles
  double total_mass = surfacedensity * r->boxsize.x * r->boxsize.y;
  double mass = 0;
  int i = 0;
  while (mass < total_mass && i < numParticles) {
    // Make all the particles the same size and touch eachother (size 2)
    struct reb_particle pt;
    if (i < numParticles / 8) {
      // spread 1/4th of the particles out along the box size
      // pt.x = (r->boxsize.x - 2 * r->boxsize.x * (i & 1)) / (i/4 + 2);
      // pt.y = (r->boxsize.y - 2 * r->boxsize.y * ((i & 2) >> 1)) / (i/4 + 2);
      // pt.z = (r->boxsize.z - 2 * r->boxsize.z * ((i & 4) >> 2)) / (i/4 + 2);
      pt.x = (r->boxsize.x / 2) / (i + 1) - 1;
      pt.y = (r->boxsize.y / 2) / (i + 1) - 1;
      pt.z = -(r->boxsize.z / 2) / (i + 1) + 1;
    } else {
      pt.x = -r->boxsize.x / 2. + i;
      pt.y = -r->boxsize.y / 2. + i;
      pt.z = 0.5; // m
    }
    pt.vx = 0;
    pt.vy = -1.5 * pt.x * OMEGA;
    pt.vz = 0;
    pt.ax = 0;
    pt.ay = 0;
    pt.az = 0;
    double radius = 2.0;
    pt.r = radius; // m
    double particle_mass =
        particle_density * 4. / 3. * M_PI * radius * radius * radius;
    pt.m = particle_mass; // kg
    reb_simulation_add(r, pt);
    mass += particle_mass;
    i += 1;
  }
  return r;
}

void mockRandomShearingSheetSimulation(struct reb_simulation **r1out,
                                       struct reb_simulation **r2out,
                                       int collisionDetectionType,
                                       int numParticles, double boxsizeIn,
                                       int ghostx, int ghosty, int ghostz) {
  struct reb_simulation *r1 = reb_simulation_create();
  struct reb_simulation *r2 = reb_simulation_create();
  // Setup constants
  r1->opening_angle2 = .5; // This determines the precission of the tree code
                           // gravity calculation.
  r1->integrator = REB_INTEGRATOR_SEI;
  r1->boundary = REB_BOUNDARY_SHEAR;
  r1->gravity = REB_GRAVITY_TREE;
  r1->collision = REB_COLLISION_TREE;
  r1->collision_resolve = reb_collision_resolve_hardsphere;
  double OMEGA = 0.00013143527; // 1/s
  r1->ri_sei.OMEGA = OMEGA;
  r1->G = 6.67428e-11;               // N / (1e-5 kg)^2 m^2
  r1->softening = 0.1;               // m
  r1->dt = 1e-3 * 2. * M_PI / OMEGA; // s
  r1->heartbeat = heartbeat;         // function pointer for heartbeat

  r2->opening_angle2 = .5; // This determines the precission of the tree code
                           // gravity calculation.
  r2->integrator = REB_INTEGRATOR_SEI;
  r2->boundary = REB_BOUNDARY_SHEAR;
  r2->gravity = REB_GRAVITY_TREE;
  r2->collision = REB_COLLISION_TREE;
  r2->collision_resolve = reb_collision_resolve_hardsphere;
  r2->ri_sei.OMEGA = OMEGA;
  r2->G = 6.67428e-11;               // N / (1e-5 kg)^2 m^2
  r2->softening = 0.1;               // m
  r2->dt = 1e-3 * 2. * M_PI / OMEGA; // s
  r2->heartbeat = heartbeat;         // function pointer for heartbeat
  // This example uses two root boxes in the x and y direction.
  // Although not necessary in this case, it allows for the parallelization
  // using MPI. See Rein & Liu for a description of what a root box is in this
  // context.
  double surfacedensity = 400;    // kg/m^2
  double particle_density = 400;  // kg/m^3
  double particle_radius_min = 1; // m
  double particle_radius_max = 4; // m
  double particle_radius_slope = -3;
  double boxsize = boxsizeIn; // m
  reb_simulation_configure_box(r1, boxsize, 2, 2, 1);
  reb_simulation_configure_box(r2, boxsize, 2, 2, 1);
  r1->N_ghost_x = ghostx;
  r1->N_ghost_y = ghosty;
  r1->N_ghost_z = ghostz;

  r2->N_ghost_x = ghostx;
  r2->N_ghost_y = ghosty;
  r2->N_ghost_z = ghostz;

  // Initial conditio
  // Use Bridges et al coefficient of restitution.
  r1->coefficient_of_restitution = coefficient_of_restitution_bridges;
  r2->coefficient_of_restitution = coefficient_of_restitution_bridges;
  // When two particles collide and the relative velocity is zero, the might
  // sink into each other in the next time step. By adding a small repulsive
  // velocity to each collision, we prevent this from happening.
  r1->minimum_collision_velocity =
      particle_radius_min * OMEGA *
      0.001; // small fraction of the shear accross a particle
  r2->minimum_collision_velocity =
      particle_radius_min * OMEGA *
      0.001; // small fraction of the shear accross a particle

  // Add all ring paricles
  for (int i = 0; i < numParticles; i++) {
    struct reb_particle pt1;
    struct reb_particle pt2;
    pt1.x = reb_random_uniform(r1, -r1->boxsize.x / 2., r1->boxsize.x / 2.);
    pt2.x = pt1.x;
    pt1.y = reb_random_uniform(r1, -r1->boxsize.y / 2., r1->boxsize.y / 2.);
    pt2.y = pt1.y;
    pt1.z = reb_random_normal(r1, 1.); // m
    pt2.z = pt1.z;
    pt1.vx = 0;
    pt2.vx = pt1.vx;
    pt1.vy = -1.5 * pt1.x * OMEGA;
    pt2.vy = pt1.vy;
    pt1.vz = 0;
    pt2.vz = pt1.vz;
    pt1.ax = 0;
    pt2.ax = pt1.ax;
    pt1.ay = 0;
    pt2.ay = pt1.ay;
    pt1.az = 0;
    pt2.az = pt1.az;
    double radius = reb_random_powerlaw(
        r1, particle_radius_min, particle_radius_max, particle_radius_slope);
    pt1.r = radius; // m
    pt2.r = pt1.r;
    double particle_mass =
        particle_density * 4. / 3. * M_PI * radius * radius * radius;
    pt1.m = particle_mass; // kg
    pt2.m = pt1.m;
    reb_simulation_add(r1, pt1);
    reb_simulation_add(r2, pt2);
  }
  *r1out = r1;
  *r2out = r2;
}

// Tests  ---------------------------------

void test_Simulation() {
  struct reb_simulation *r = mockShearingSheetSimulation(REB_COLLISION_TREE, 2);
  reb_simulation_integrate(r, 10);
  assert(r->collisions_N == 110);
  // printf("\n\n%d\n", r->collisions_N);
}

void test_CollisionSearchOriginal() {
  struct reb_simulation *r = mockShearingSheetSimulation(REB_COLLISION_TREE, 8);
  mock_reb_collision_search(r);
  // assert(r->collisions_N == 34);
  printf("\n\n%d\n", r->collisions_N);
}

void test_CollisionSearchCompare(int particleCount) {
  struct reb_simulation *r1 =
      mockShearingSheetSimulation(REB_COLLISION_TREE, particleCount);
  struct reb_simulation *r2 =
      mockShearingSheetSimulation(REB_COLLISION_TREE, particleCount);
  mock_reb_collision_search(r1);
  int collisions_original = r1->collisions_N;
  speedup_reb_collision_search(r2);
  int collisions_new = r2->collisions_N;
  assert(collisions_original == collisions_new);
  reb_simulation_free(r1);
  reb_simulation_free(r2);
}

void test_CollisionSearchCompareRandom(int particleCount, double boxSize,
                                       int ghostx, int ghosty, int ghostz) {
  struct reb_simulation *r1Ptr = NULL;
  struct reb_simulation *r2Ptr = NULL;
  mockRandomShearingSheetSimulation(&r1Ptr, &r2Ptr, REB_COLLISION_TREE,
                                    particleCount, boxSize, ghostx, ghosty,
                                    ghostz);
  struct reb_simulation *r1 = r1Ptr;
  struct reb_simulation *r2 = r2Ptr;
  mock_reb_collision_search(r1);
  int collisions_original = r1->collisions_N;
  speedup_reb_collision_search(r2);
  int collisions_new = r2->collisions_N;
  assert(collisions_original == collisions_new);
  reb_simulation_free(r1);
  reb_simulation_free(r2);
}

void time_CollisionSearchOriginal() {
  double total_time = 0.0;
  struct reb_simulation *r;
  for (int i = 0; i < 50000; i++) {
    r = mockShearingSheetSimulation(REB_COLLISION_TREE, i % 200);
    total_time += mock_reb_collision_search(r);
    reb_simulation_free(r);
  }
  // assert(r->collisions_N == 110);
  printf("\nTotal Time: %f\n", total_time);
}

void time_CollisionSearchSpeedup() {
  double total_time = 0.0;
  struct reb_simulation *r;
  for (int i = 0; i < 50000; i++) {
    r = mockShearingSheetSimulation(REB_COLLISION_TREE, i % 200);
    total_time += speedup_reb_collision_search(r);
    reb_simulation_free(r);
  }
  printf("\nTotal Time: %f\n", total_time);
}

void time_CollisionSearchSpeedupCompare() {
  double total_time1 = 0.0;
  double total_time2 = 0.0;
  struct reb_simulation *r;
  for (int i = 0; i < 50000; i++) {
    r = mockShearingSheetSimulation(REB_COLLISION_TREE, i % 200);
    total_time1 += mock_reb_collision_search(r);
    reb_simulation_free(r);
  }

  for (int i = 0; i < 50000; i++) {
    r = mockShearingSheetSimulation(REB_COLLISION_TREE, i % 200);
    total_time2 += speedup_reb_collision_search(r);
    reb_simulation_free(r);
  }
  printf("\nTotal Time Original: %f\n", total_time1);
  printf("Total Time Speedup: %f\n", total_time2);
  printf("Relative Time %f\n", (total_time2 / total_time1));
}

typedef struct timeResults {
  double timeOrig;
  double timeNew;
} timeResults;

void time_CollisionSearchCompareRandom(int particleCount, double boxSize,
                                       int ghostx, int ghosty, int ghostz,
                                       int iterations,
                                       struct timeResults *timeResults) {
  double total_time1 = 0.0;
  double total_time2 = 0.0;
  printf("Testing with ParticleCount: %d, BoxSize: %f, ghostx: %d, ghostz: %d, "
         "ghostz: %d, iterations: %d\n",
         particleCount, boxSize, ghostx, ghosty, ghostz, iterations);
  for (int i = 0; i < iterations; i++) {
    struct reb_simulation *r1Ptr = NULL;
    struct reb_simulation *r2Ptr = NULL;
    mockRandomShearingSheetSimulation(&r1Ptr, &r2Ptr, REB_COLLISION_TREE,
                                      particleCount, boxSize, ghostx, ghosty,
                                      ghostz);
    struct reb_simulation *r1 = r1Ptr;
    struct reb_simulation *r2 = r2Ptr;
    total_time1 += mock_reb_collision_search(r1);
    total_time2 += speedup_reb_collision_search(r2);
    mock_reb_collision_search(r1);
    speedup_reb_collision_search(r2);
    reb_simulation_free(r1);
    reb_simulation_free(r2);
  }
  printf("\nTotal Time Original: %f\n", total_time1);
  printf("Total Time Speedup: %f\n", total_time2);
  printf("Relative Time %f\n", (total_time2 / total_time1));

  timeResults->timeOrig += total_time1;
  timeResults->timeNew += total_time2;
}

int main(int argc, char *argv[]) {

  // Choose collision method (Case section in collision.c)
  // mockShearingSheetSimulation(Collision type (Keep REB_COLLISION_TREE,
  // numParticles)) reb_simulation_start_server();
  //  mock_tree_collisions();

  // TEST
  test_CollisionSearchCompare(-1);
  test_CollisionSearchCompare(0);
  test_CollisionSearchCompare(1);
  test_CollisionSearchCompare(2);
  test_CollisionSearchCompare(10);
  test_CollisionSearchCompare(100);
  int numParticles[] = {
      0, 1, 10, 100, 1000,
  };
  double boxSize[] = {0.0, 1.0, 10.0, 100.0, 1000.0};
  int ghostx[] = {0, 1, 2, 5, 10};
  int ghosty[] = {0, 1, 2, 5, 10};
  int ghostz[] = {0, 1, 2, 5, 10};

  printf(
      "Running 3125 Test for \nNumber of Particles: 0,1,10,100,1000\nBoxSize: "
      "0.0, 1.0, 10.0, 100.0, 1000.0\nghost x, z, and y: 0,1,2,5,10\n");
  for (int a = 0; a < 5; a++) {
    for (int b = 0; b < 5; b++) {
      for (int c = 0; c < 5; c++) {
        for (int d = 0; d < 5; d++) {
          for (int e = 0; e < 5; e++) {
            test_CollisionSearchCompareRandom(numParticles[a], boxSize[b],
                                              ghostx[c], ghosty[d], ghostz[e]);
          }
        }
      }
    }
  }
  // test_CollisionSearchCompareRandom(100, 100.0, 2, 2, 0);

  // test_CollisionSearchCompare(1000);
  // WANT TO ADD END TO END TEST HERE WITH OUTPUT TO BIN

  // TIME
  //  time_CollisionSearchSpeedup();
  // time_CollisionSearchSpeedupCompare();
  struct timeResults *timeresults = malloc(sizeof(struct timeResults));
  printf(
      "Running 3125 Test for \nNumber of Particles: 0,1,10,100,1000\nBoxSize: "
      "0.0, 1.0, 10.0, 100.0, 1000.0\nghost x, z, and y: 0,1,2,5,10\n");
  int numParticles2[] = {100, 1000};
  double boxSize2[] = {10.0, 100.0};
  int ghostx2[] = {1, 2};
  int ghosty2[] = {1, 2};
  int ghostz2[] = {1, 2};
  int iterations = 15;
  for (int a = 0; a < 2; a++) {
    for (int b = 0; b < 2; b++) {
      for (int c = 0; c < 2; c++) {
        for (int d = 0; d < 2; d++) {
          for (int e = 0; e < 2; e++) {
            time_CollisionSearchCompareRandom(
                numParticles2[a], boxSize2[b], ghostx2[c], ghosty2[d],
                ghostz2[e], iterations, timeresults);
          }
        }
      }
    }
  }
  printf("================================\n");
  printf("Total Relative Time %f\n",
         (timeresults->timeNew / timeresults->timeOrig));
  printf("================================\n");

  // test_Simulation(r, 10);

  // struct reb_simulation* r = mockShearingSheetSimulation(REB_COLLISION_TREE,
  // 8); time_FullCollisionSearchOriginal(r);
  // time_CollisionSearchOriginal(r);
  // reb_simulation_free(r);
  free(timeresults);
}
