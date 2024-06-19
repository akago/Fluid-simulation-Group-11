#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include "Solver.h"
#include "util.h"
#include "ConstraintSolver.h"


/* macro */
#define IX(i,j) ((i)+(N+2)*(j))
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
#define END_FOR }}


/*
	Particle Sovler
*/
void UpdateForces(std::vector<Particle*> pVector, std::vector<RigidBody*> rigid_bodies, std::vector<Force*> fVector) {
	int ii, size = pVector.size();

	// Clear forces
	for (ii = 0; ii < size; ii++) {
		pVector[ii]->clearForce();
	}

	for (auto rb : rigid_bodies) {
		rb->clearForce();
	}

	// Compute forces
	for (auto force : fVector) {
		force->applyForce();
	}

	// Update rigid bodies forces
	for (auto rb : rigid_bodies) {
		rb->addForce();
		rb->addTorque();
	}
}



void ParticleDeriv(std::vector<Particle*> pVector, std::vector<RigidBody*> rigid_bodies, std::vector<Force*> fVector, std::vector<Constraint*> cVector, double* dst, double dt) {
	UpdateForces(pVector, rigid_bodies, fVector);
	ApplyConstraintForce(pVector, cVector);

	for (auto p : pVector) {
		*(dst++) = p->m_Velocity[0];
		*(dst++) = p->m_Velocity[1];
		*(dst++) = p->m_Force[0] / p->mass;
		*(dst++) = p->m_Force[1] / p->mass;
	}

	/* rigid body simulation*/
	for (int i = 0; i < rigid_bodies.size(); i++) {
		rigid_bodies[i]->m_Velocity += (rigid_bodies[i]->m_Force / rigid_bodies[i]->Mass) * dt;
		rigid_bodies[i]->m_Angular_Vec += (rigid_bodies[i]->m_torque / rigid_bodies[i]->m_Inertia) * dt;
		rigid_bodies[i]->m_Position += rigid_bodies[i]->m_Velocity * dt;
		rigid_bodies[i]->m_Angle += rigid_bodies[i]->m_Angular_Vec * dt;

		rigid_bodies[i]->updateR(rigid_bodies[i]->m_Angle);

		rigid_bodies[i]->m_Angular_Vec *= 0.99f;

		/* Update vertices local position*/
		for (int j = 0; j < rigid_bodies[i]->m_Vertices.size(); j++)
		{
			rigid_bodies[i]->m_Vertices[j]->m_Position = rigid_bodies[i]->m_R.multiByVec2(rigid_bodies[i]->m_Vertices[j]->m_Position);
		}

		
	}

}

void GetSystemState(std::vector<Particle*> pVector, double* dst) {
	for (auto p : pVector) {
		*(dst++) = p->m_Position[0];
		*(dst++) = p->m_Position[1];
		*(dst++) = p->m_Velocity[0];
		*(dst++) = p->m_Velocity[1];
	}
}

void SetSystemState(std::vector<Particle*> pVector, double* src) {
	for (auto p : pVector) {
		p->m_Position[0] = *(src++);
		p->m_Position[1] = *(src++);
		p->m_Velocity[0] = *(src++);
		p->m_Velocity[1] = *(src++);
	}
}

double* ComputeK1(std::vector<Particle*> pVector, std::vector<RigidBody*> rigid_bodies, std::vector<Force*> fVector, std::vector<Constraint*> cVector, double dt) {
	double* k1 = (double*)malloc(pVector.size() * 4 * sizeof(double));

	// Compute k1
	ParticleDeriv(pVector, rigid_bodies, fVector, cVector, k1, dt);
	vecTimesScalar(pVector.size() * 4, k1, dt);

	return k1;
}

double* ComputeKHelper(std::vector<Particle*> pVector, std::vector<RigidBody*> rigid_bodies, std::vector<Force*> fVector, std::vector<Constraint*> cVector, double dt, double* k, double factor) {
	double* k_next = (double*)malloc(pVector.size() * 4 * sizeof(double));

	// Compute next k
	GetSystemState(pVector, k_next);
	vecAddEqualWithFactor(pVector.size() * 4, k_next, k, factor);

	SetSystemState(pVector, k_next);
	ParticleDeriv(pVector, rigid_bodies, fVector, cVector, k_next, dt);
	vecTimesScalar(pVector.size() * 4, k_next, dt);

	return k_next;

}

double* ComputeK2(std::vector<Particle*> pVector, std::vector<RigidBody*> rigid_bodies, std::vector<Force*> fVector, std::vector<Constraint*> cVector, double dt, double* k1) {
	return ComputeKHelper(pVector, rigid_bodies, fVector, cVector, dt, k1, 0.5);
}

double* ComputeK3(std::vector<Particle*> pVector, std::vector<RigidBody*> rigid_bodies, std::vector<Force*> fVector, std::vector<Constraint*> cVector, double dt, double* k2) {
	return ComputeKHelper(pVector, rigid_bodies, fVector, cVector, dt, k2, 0.5);
}

double* ComputeK4(std::vector<Particle*> pVector, std::vector<RigidBody*> rigid_bodies, std::vector<Force*> fVector, std::vector<Constraint*> cVector, double dt, double* k3) {
	return ComputeKHelper(pVector, rigid_bodies, fVector, cVector, dt, k3, 1);
}

void Euler_step(std::vector<Particle*> pVector, std::vector<RigidBody*> rigid_bodies, std::vector<Force*> fVector, std::vector<Constraint*> cVector, float dt) {
	double* state = (double*)malloc(pVector.size() * 4 * sizeof(double));
	GetSystemState(pVector, state);

	double* k1 = ComputeK1(pVector, rigid_bodies, fVector, cVector, dt);
	vecAddEqual(pVector.size() * 4, state, k1);
	SetSystemState(pVector, state);

	free(state);
	free(k1);
}

void Midpoint_step(std::vector<Particle*> pVector, std::vector<RigidBody*> rigid_bodies, std::vector<Force*> fVector, std::vector<Constraint*> cVector, float dt) {
	double* state = (double*)malloc(pVector.size() * 4 * sizeof(double));
	GetSystemState(pVector, state);

	double* k1 = ComputeK1(pVector, rigid_bodies, fVector, cVector, dt);
	double* k2 = ComputeK2(pVector, rigid_bodies, fVector, cVector, dt, k1);

	vecAddEqual(pVector.size() * 4, state, k2);
	SetSystemState(pVector, state);

	free(state);
	free(k1);
	free(k2);
}

void Runge_Kutta_4(std::vector<Particle*> pVector, std::vector<RigidBody*> rigid_bodies, std::vector<Force*> fVector, std::vector<Constraint*> cVector, float dt) {
	double* state = (double*)malloc(pVector.size() * 4 * sizeof(double));
	GetSystemState(pVector, state);

	double* k1 = ComputeK1(pVector, rigid_bodies, fVector, cVector, dt);
	double* k2 = ComputeK2(pVector, rigid_bodies, fVector, cVector, dt, k1);
	SetSystemState(pVector, state);
	double* k3 = ComputeK3(pVector, rigid_bodies, fVector, cVector, dt, k2);
	SetSystemState(pVector, state);
	double* k4 = ComputeK4(pVector, rigid_bodies, fVector, cVector, dt, k3);

	vecAddEqualWithFactor(pVector.size() * 4, state, k1, 1.0 / 6);
	vecAddEqualWithFactor(pVector.size() * 4, state, k2, 1.0 / 3);
	vecAddEqualWithFactor(pVector.size() * 4, state, k3, 1.0 / 3);
	vecAddEqualWithFactor(pVector.size() * 4, state, k4, 1.0 / 6);
	SetSystemState(pVector, state);

	free(state);
	free(k1);
	free(k2);
	free(k3);
	free(k4);
}



/*
	Fluid Sovler
*/

void add_source ( int N, float * x, float * s, float dt )
{
	int i, size=(N+2)*(N+2);
	for ( i=0 ; i<size ; i++ ) x[i] += dt*s[i];
}


void set_intbnd(int N, float * x, int* mask, int i, int j, int b){
	if (!mask[IX(i+1, j)]) x[IX(i,j)] = b==1 ? -x[IX(i+1, j)] : x[IX(i+1, j)];
	else if (!mask[IX(i-1, j)]) x[IX(i,j)] = b==1 ? -x[IX(i-1, j)] : x[IX(i-1, j)];
	else if (!mask[IX(i, j+1)]) x[IX(i,j)] = b==2 ? -x[IX(i, j+1)] : x[IX(i, j+1)];
	else if (!mask[IX(i, j-1)]) x[IX(i,j)] = b==2 ? -x[IX(i, j-1)] : x[IX(i, j-1)];
	else if (!mask[IX(i-1, j-1)]) x[IX(i,j)] = 0.5f*(x[IX(i-1,j)]+x[IX(i,j-1)]);
	else if (!mask[IX(i+1, j-1)]) x[IX(i,j)] = 0.5f*(x[IX(i+1,j)]+x[IX(i,j-1)]);
	else if (!mask[IX(i-1, j+1)]) x[IX(i,j)] = 0.5f*(x[IX(i-1,j)]+x[IX(i,j+1)]);
	else if (!mask[IX(i+1, j+1)]) x[IX(i,j)] = 0.5f*(x[IX(i+1,j)]+x[IX(i,j+1)]);
}

void set_bnd_vel(int N, float * x, int* mask, float* bnd_vel, int i, int j, int b) {
	if (!mask[IX(i + 1, j)]) x[IX(i, j)] = b == 1 ? bnd_vel[IX(i, j)]-x[IX(i + 1, j)] : x[IX(i + 1, j)];
	else if (!mask[IX(i - 1, j)]) x[IX(i, j)] = b == 1 ? bnd_vel[IX(i, j)] -x[IX(i - 1, j)] : x[IX(i - 1, j)];
	else if (!mask[IX(i, j + 1)]) x[IX(i, j)] = b == 2 ? bnd_vel[IX(i, j)] -x[IX(i, j + 1)] : x[IX(i, j + 1)];
	else if (!mask[IX(i, j - 1)]) x[IX(i, j)] = b == 2 ? bnd_vel[IX(i, j)] -x[IX(i, j - 1)] : x[IX(i, j - 1)];
	else if (!mask[IX(i - 1, j - 1)]) x[IX(i, j)] = 0.5f*(x[IX(i - 1, j)] + x[IX(i, j - 1)]);
	else if (!mask[IX(i + 1, j - 1)]) x[IX(i, j)] = 0.5f*(x[IX(i + 1, j)] + x[IX(i, j - 1)]);
	else if (!mask[IX(i - 1, j + 1)]) x[IX(i, j)] = 0.5f*(x[IX(i - 1, j)] + x[IX(i, j + 1)]);
	else if (!mask[IX(i + 1, j + 1)]) x[IX(i, j)] = 0.5f*(x[IX(i + 1, j)] + x[IX(i, j + 1)]);
}

void set_bnd ( int N, int b, float * x, int* mask, float* bnd_vel)
{
	int i,j;

	for ( i=1 ; i<=N ; i++ ) {
		x[IX(0  ,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)];
		x[IX(N+1,i)] = b==1 ? -x[IX(N,i)] : x[IX(N,i)];
		x[IX(i,0  )] = b==2 ? -x[IX(i,1)] : x[IX(i,1)];
		x[IX(i,N+1)] = b==2 ? -x[IX(i,N)] : x[IX(i,N)];
	}
	x[IX(0  ,0  )] = 0.5f*(x[IX(1,0  )]+x[IX(0  ,1)]);
	x[IX(0  ,N+1)] = 0.5f*(x[IX(1,N+1)]+x[IX(0  ,N)]);
	x[IX(N+1,0  )] = 0.5f*(x[IX(N,0  )]+x[IX(N+1,1)]);
	x[IX(N+1,N+1)] = 0.5f*(x[IX(N,N+1)]+x[IX(N+1,N)]);

	/* internal boundary*/
	FOR_EACH_CELL
		if (mask[IX(i, j)] != 0){
			if (b == 0) {
				set_intbnd(N, x, mask, i, j, b);
			}
			else {
				set_bnd_vel(N, x, mask, bnd_vel, i, j, b);
			}
		}
		
	END_FOR
	FOR_EACH_CELL
		if (mask[IX(i, j)] != 0) {
			if (b == 0) {
				set_intbnd(N, x, mask, i, j, b);
			}
			else {
				set_bnd_vel(N, x, mask, bnd_vel, i, j, b);
			}

		}
	END_FOR

}

void lin_solve ( int N, int b, float * x, float * x0, int* mask, float* bnd_vel, float a, float c )
{
	int i, j, k;

	for ( k=0 ; k<20 ; k++ ) {
		FOR_EACH_CELL
			x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/c;
		END_FOR
		set_bnd ( N, b, x, mask, bnd_vel);
	}
}

void diffuse ( int N, int b, float * x, float * x0, int* mask, float* bnd_vel, float diff, float dt )
{
	float a=dt*diff*N*N; // 1/h = N
	lin_solve ( N, b, x, x0, mask, bnd_vel, a, 1+4*a );
}

void advect ( int N, int b, float * d, float * d0, float * u, float * v, int* mask, float* bnd_vel, float dt )
{
	int i, j, i0, j0, i1, j1;
	float x, y, s0, t0, s1, t1, dt0;

	dt0 = dt*N;
	FOR_EACH_CELL
		x = i-dt0*u[IX(i,j)]; y = j-dt0*v[IX(i,j)];
		if (x<0.5f) x=0.5f; if (x>N+0.5f) x=N+0.5f; i0=(int)x; i1=i0+1;
		if (y<0.5f) y=0.5f; if (y>N+0.5f) y=N+0.5f; j0=(int)y; j1=j0+1;
		s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
		d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+
					 s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
	END_FOR
	set_bnd ( N, b, d, mask, bnd_vel);
}

void project ( int N, float * u, float * v, float * p, float * div, int* mask, float* bnd_vel_u,  float* bnd_vel_v)
{
	int i, j;

	FOR_EACH_CELL
		div[IX(i,j)] = -0.5f*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)])/N;
		p[IX(i,j)] = 0;
	END_FOR	
	set_bnd ( N, 0, div, mask, bnd_vel_u); set_bnd ( N, 0, p, mask, bnd_vel_u);

	lin_solve ( N, 0, p, div, mask, bnd_vel_u, 1, 4 );

	FOR_EACH_CELL
		u[IX(i,j)] -= 0.5f*N*(p[IX(i+1,j)]-p[IX(i-1,j)]);
		v[IX(i,j)] -= 0.5f*N*(p[IX(i,j+1)]-p[IX(i,j-1)]);
	END_FOR
	set_bnd ( N, 1, u, mask, bnd_vel_u); set_bnd ( N, 2, v, mask, bnd_vel_v);
}



void vorticity_confinement(int N, float * u, float * v, float * u0, float * v0, float epsilon){
	int i, j, k;
	float* omega, *eta_x, *eta_y;
	int size = (N+2)*(N+2);
	
	omega = (float *) malloc ( size*sizeof(float) );
	eta_x = (float *) malloc ( size*sizeof(float) );
	eta_y = (float *) malloc ( size*sizeof(float) );


	//ω= ∂v/∂x - ∂u/∂y
	FOR_EACH_CELL
			omega[IX(i,j)] = 0.5f*N*(v[IX(i+1,j)]-v[IX(i-1,j)] - u[IX(i,j+1)]-u[IX(i,j-1)]);
	END_FOR

	// ∇∣ω∣=(∂∣ω∣/∂x, ∂∣ω∣/∂y)
	FOR_EACH_CELL
			eta_x[IX(i,j)] = 0.5f*N* (fabs(omega[IX(i+1,j)]) - fabs(omega[IX(i-1,j)]));
            eta_y[IX(i,j)] = 0.5f*N* (fabs(omega[IX(i,j+1)]) - fabs(omega[IX(i,j-1)]));
	END_FOR

	// N = eta / |eta|
	FOR_EACH_CELL
			float magnitude = sqrt(eta_x[IX(i,j)] * eta_x[IX(i,j)] + eta_y[IX(i,j)] * eta_y[IX(i,j)]);
			if (magnitude > 0) {
                eta_x[IX(i,j)] /= magnitude;
				eta_y[IX(i,j)] /= magnitude;

				u0[IX(i,j)] += epsilon * (eta_y[IX(i,j)] * omega[IX(i,j)]) / N;		
				v0[IX(i,j)] += -epsilon * (eta_x[IX(i,j)] * omega[IX(i,j)]) / N;		
            }
	END_FOR

	free(omega);
	free(eta_x);
	free(eta_y);

}


void dens_step ( int N, float * x, float * x0, float * u, float * v, int * mask, float* bnd_vel_u,  float* bnd_vel_v, float diff, float dt )
{
	add_source ( N, x, x0, dt );
	SWAP ( x0, x ); diffuse ( N, 0, x, x0, mask, bnd_vel_u, diff, dt );
	SWAP ( x0, x ); advect ( N, 0, x, x0, u, v, mask, bnd_vel_u, dt );
}

void vel_step ( int N, float * u, float * v, float * u0, float * v0, int * mask, float* bnd_vel_u,  float* bnd_vel_v, float visc, float dt )
{	
	vorticity_confinement(N, u, v, u0, v0, 0.1);
	add_source ( N, u, u0, dt ); add_source ( N, v, v0, dt );
	SWAP ( u0, u ); diffuse ( N, 1, u, u0, mask, bnd_vel_u, visc, dt );
	SWAP ( v0, v ); diffuse ( N, 2, v, v0, mask, bnd_vel_v, visc, dt );
	project ( N, u, v, u0, v0, mask, bnd_vel_u,  bnd_vel_v);
	SWAP ( u0, u ); SWAP ( v0, v );
	advect ( N, 1, u, u0, u0, v0, mask, bnd_vel_u, dt ); advect ( N, 2, v, v0, u0, v0, mask, bnd_vel_v, dt );
	project ( N, u, v, u0, v0, mask, bnd_vel_u,  bnd_vel_v);
}




