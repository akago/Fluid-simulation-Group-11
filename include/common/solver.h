#include <stdlib.h>
#include <stdbool.h>
#include <vector>

#include "RigidBody.h"
#include "Particle.h"
#include "Force.h"
#include "Constraint.h"
#include "linearSolver.h"

#define DAMP 0.98f
#define RAND (((rand()%2000)/1000.f)-1.f)

/*
	Particle solver
*/

enum IntegrationType
{
	Euler,
	Midpoint,
	RungeKutta
};

typedef void(*IntegrationFunctionHook)(std::vector<Particle*>, std::vector<RigidBody*>, std::vector<Force*>, std::vector<Constraint*>, float);

void ParticleDeriv(std::vector<Particle*> pVector, std::vector<RigidBody*> rigid_bodies, std::vector<Force*> fVector, std::vector<Constraint*> cVector, double* dst, double dt);
void GetSystemState(std::vector<Particle*> pVector, double* dst);
void SetSystemState(std::vector<Particle*> pVector, double* src);

void Euler_step(std::vector<Particle*> pVector, std::vector<RigidBody*> rigid_bodies, std::vector<Force*> fVector, std::vector<Constraint*> cVector, float dt);
void Midpoint_step(std::vector<Particle*> pVector, std::vector<RigidBody*> rigid_bodies, std::vector<Force*> fVector, std::vector<Constraint*> cVector, float dt);
void Runge_Kutta_4(std::vector<Particle*> pVector, std::vector<RigidBody*> rigid_bodies, std::vector<Force*> fVector, std::vector<Constraint*> cVector, float dt);


/*
	Fluid solver
*/


void add_source ( int N, float * x, float * s, float dt );
void set_bnd (  int N, int b, float * x, int* mask, float* bnd_vel);
void lin_solve ( int N, int b, float * x, float * x0, int* mask, float* bnd_vel, float a, float c  );
void diffuse ( int N, int b, float * x, float * x0, int* mask, float* bnd_vel, float diff, float dt );
void advect ( int N, int b, float * d, float * d0, float * u, float * v, int* mask, float* bnd_vel, float dt );
void project ( int N, float * u, float * v, float * p, float * div, int* mask, float* bnd_vel_u,  float* bnd_vel_v);
void vorticity_confinement(int N, float * u, float * v, float * u0, float * v0, float epsilon);


void dens_step ( int N, float * x, float * x0, float * u, float * v, int * mask, float* bnd_vel_u,  float* bnd_vel_v, float diff, float dt );
void vel_step ( int N, float * u, float * v, float * u0, float * v0, int * mask, float* bnd_vel_u,  float* bnd_vel_v, float visc, float dt );
