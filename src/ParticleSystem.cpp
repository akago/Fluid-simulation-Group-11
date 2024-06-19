#include "ParticleSystem.h"

ParticleSystem::ParticleSystem() {
    integrationMethod = Euler_step;
}

/*
----------------------------------------------------------------------
Set hook function -- plugable integration scheme
Euler: Euler step
Midpoint: Midpoint step
RungeKutta: RungeKutta-4 step
----------------------------------------------------------------------
*/
void ParticleSystem::setIntegrationHook(IntegrationType t) {
    switch (t)
	{
		case Euler:
			integrationMethod = Euler_step;
			break;
		case Midpoint:
			integrationMethod = Midpoint_step;
			break;
		case RungeKutta:
			integrationMethod = Runge_Kutta_4;
			break;
		default:
			integrationMethod = Euler_step; //default approach
			break;
	}
}

void ParticleSystem::setDt(float dt) {
    this->dt = dt;
}

void ParticleSystem::simulationStep() {
	for(int i = 0; i < particles.size(); i++) {
		previousPositions[i] = particles[i]->m_Position;
	}

    integrationMethod(particles, rigid_bodies, forces, constraints, dt);
	
	for(int i = 0; i < particles.size(); i++) {
		for(Wall* wall : walls) {
			wall->adjustParticle(particles[i], previousPositions[i]);
		}
	}
}

void ParticleSystem::projectRigidBodies(int N, int* internal_bd, float* bnd_vel_u, float* bnd_vel_v) {
	for (int i = 0; i < rigid_bodies.size(); i++) {
		rigid_bodies[i]->projectToGrid(N, internal_bd, bnd_vel_u, bnd_vel_v);
	}
}


void ParticleSystem::addParticle(Particle* particle) {
    particles.push_back(particle);
	previousPositions.push_back(Vec2f());
}

void ParticleSystem::addForce(Force* force) {
    forces.push_back(force);
}

void ParticleSystem::addConstraint(Constraint* constraint) {
    constraints.push_back(constraint);
}

void ParticleSystem::addWall(Wall* wall) {
	walls.push_back(wall);
}

void ParticleSystem::addRigid(RigidBody* rigid) {
	rigid_bodies.push_back(rigid);
}


void ParticleSystem::removeLastForce() {
	forces.pop_back();
}

std::vector<Particle*>& ParticleSystem::getParticles() {
    return particles;
}
std::vector<Force*>& ParticleSystem::getForces() {
    return forces;
}
std::vector<Constraint*>& ParticleSystem::getConstraints() {
    return constraints;
}
std::vector<Wall*>& ParticleSystem::getWalls() {
	return walls;
}

std::vector<RigidBody*>& ParticleSystem::getRigids() {
	return rigid_bodies;
}

int ParticleSystem::particleCount() {
    return particles.size();
}

void ParticleSystem::drawParticles ( )
{
    for (auto particle : particles) {
		particle->draw();
	}
}

void ParticleSystem::drawForces ( )
{
	for (auto force : forces) {
		force->draw();
	}	
}

void ParticleSystem::drawConstraints ( )
{
	for (auto constraint : constraints) {
		constraint->draw();
	}
}

void ParticleSystem::drawWalls ( )
{
	for (auto wall : walls) {
		wall->draw();
	}		
}

void ParticleSystem::drawRigids() {
	for (auto rb : rigid_bodies) {
		rb->draw();
	}
}


void ParticleSystem::reset() {
    for (auto particle : particles) {
		particle->reset();
	}
	for (auto rb : rigid_bodies) {
		rb->reset();
	}
}

ParticleSystem::~ParticleSystem() {
	walls.clear();
	previousPositions.clear();
    particles.clear();
    forces.clear();
    constraints.clear();
}