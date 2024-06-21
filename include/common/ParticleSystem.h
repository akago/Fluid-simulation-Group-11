#pragma once

#include "Particle.h"
#include "Force.h"
#include "Solver.h"
#include "Constraint.h"
#include "RigidBody.h"
#include "Rectangle.h"
#include "Wall.h"

#include <vector>

class ParticleSystem {
    public:
        ParticleSystem();
        virtual ~ParticleSystem();
        
        void setIntegrationHook(IntegrationType t);
        virtual void simulationStep();
		void projectRigidBodies(int N, int* internal_bd, float* bnd_vel_u, float* bnd_vel_v);

        void setDt(float dt);

        void addParticle(Particle* particle);
        void addForce(Force* force);
        void addConstraint(Constraint* constraint);
		void addRigid(RigidBody* rigid);
        void addWall(Wall* wall);

        void removeLastForce();

        std::vector<Particle*>& getParticles();
        std::vector<Force*>& getForces();
        std::vector<Constraint*>& getConstraints();
        std::vector<Wall*>& getWalls();
		std::vector<RigidBody*>& getRigids();

        int particleCount();
        int objectCount();
        
        void drawParticles();
        void drawForces();
        void drawConstraints();
        void drawWalls();
		void drawRigids();

        void reset();
    protected:
        float dt = 0.1f;

        std::vector<Vec2f> previousPositions;
        std::vector<Wall*> walls;
        std::vector<Particle*> particles;
        std::vector<Force*> forces;
        std::vector<Constraint*> constraints;
		std::vector<RigidBody*> rigid_bodies;

        IntegrationFunctionHook integrationMethod;
};