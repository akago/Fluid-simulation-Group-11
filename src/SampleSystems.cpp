#include "SampleSystems.h"
#include "ParticleSystem.h"
#include "ClothParticleSystem.h"

#define PI 3.1415926535897932384626433832795

/*
	Circle, Rod, Linear Spring
*/
ParticleSystem* system1() {
    ParticleSystem* system = new ParticleSystem();

    system->setDt(0.005);

    const double dist = 0.2;
	const Vec2f center(0.0, 0.0);
	const Vec2f offset(dist, 0.0);

	// Set integration scheme.
	system->setIntegrationHook(Midpoint);

	// Create three particles, attach them to each other, then add a
	// circular wire constraint to the first.

    Particle* particle0 = new Particle(center + offset);
    Particle* particle1 = new Particle(center + offset + offset);
    Particle* particle2 = new Particle(center + offset + offset + offset);

	system->addParticle(particle0);
	system->addParticle(particle1);
	system->addParticle(particle2);
	
	system->addForce(new GravityForce(system->getParticles(), 0.01));
	system->addForce(new SpringForce(particle0, particle1, dist, 0.7, 0.15));

	// Create Global Constraint Jacobian Matrix
	Constraint::GlobalJ = new GlobalMatrix(0,system->particleCount()*2);
	Constraint::GlobalJdot = new GlobalMatrix(0, system->particleCount() * 2);
	Constraint::global_cons_num = 0;
	Constraint::kd = 0.2;
	Constraint::ks = 0.3;

	system->addConstraint(new CircularWireConstraint(0, particle0, center, dist));
	system->addConstraint(new RodConstraint(1, 2, particle1, particle2, dist));

    return system;
}

/*
	Line constraint
*/
ParticleSystem* system2() {
    ParticleSystem* system = new ParticleSystem();

    system->setDt(0.02);

	const Vec2f center(0.0, 0.0);

	// Set integration scheme.
	system->setIntegrationHook(RungeKutta);

    Particle* particle0 = new Particle(center);

	system->addParticle(particle0);
	
	system->addForce(new GravityForce(system->getParticles(), 0.005));

	// Create Global Constraint Jacobian Matrix
	Constraint::GlobalJ = new GlobalMatrix(0,system->particleCount()*2);
	Constraint::GlobalJdot = new GlobalMatrix(0, system->particleCount() * 2);
	Constraint::global_cons_num = 0;
	Constraint::kd = 0.2;
	Constraint::ks = 0.3;

	system->addConstraint(new LineWireConstraint(0, particle0, 1, -1, 0));

    return system;
}

/*
	Angular string
*/
ParticleSystem* system3() {
	ParticleSystem* system = new ParticleSystem();

	const double dist = 0.2;
	const Vec2f center(0.0, 0.0);
	const Vec2f offset(dist, 0.0);
	double alpha = degreesToRadians(120); // degrees

	// Set integration scheme.
	system->setIntegrationHook(Euler);
	system->setDt(0.02);
	// Create three particles, attach them to each other

	Particle* particle0 = new Particle(center + offset);
    Particle* particle1 = new Particle(center + offset + offset);
    Particle* particle2 = new Particle(Vec2f(center[0] + dist + dist + dist / 2, center[0] - sqrt(3) / 2 * dist));
	
	system->addParticle(particle0);
	system->addParticle(particle1);
	system->addParticle(particle2);

	system->addForce(new GravityForce(system->getParticles(), 0.001));
	system->addForce(new SpringForce(particle0, particle1, dist, 0.07, 0.15));
	system->addForce(new SpringForce(particle1, particle2, dist, 0.1, 0.15));
	system->addForce(new AngularSpring(particle0, particle1, particle2, alpha, 0.05, 0.2));


	// Create Global Constraint Jacobian Matrix
	Constraint::GlobalJ = new GlobalMatrix(0,system->particleCount()*2);
	Constraint::GlobalJdot = new GlobalMatrix(0, system->particleCount() * 2);
	Constraint::global_cons_num = 0;
	Constraint::kd = 0.2;
	Constraint::ks = 0.3;

	system->addConstraint(new CircularWireConstraint(0, particle0, center, dist));

	return system;
}


ParticleSystem* system4() {
	ParticleSystem* system = new ParticleSystem();

	const Vec2f center(0.0, 0.0);

	// Set integration scheme.
	system->setIntegrationHook(Euler);
	system->setDt(0.02);

	// Create three particles, attach them to each other

	Particle* particle0 = new Particle(center);
	
	system->addParticle(particle0);

	system->addForce(new GravityForce(system->getParticles(), 0.05));

	system->addWall(new Wall(Vec2f(0,-0.5), 0.1, 0.4));

	// Create Global Constraint Jacobian Matrix
	Constraint::GlobalJ = new GlobalMatrix(0,system->particleCount()*2);
	Constraint::GlobalJdot = new GlobalMatrix(0, system->particleCount() * 2);
	Constraint::global_cons_num = 0;
	Constraint::kd = 0.2;
	Constraint::ks = 0.3;

	return system;
}

ParticleSystem* cloth1() {
	double posX = -0.5;
    double posY = 0.5;
    int width = 10;
    int height = 10;
    double dist = 0.1;
    double mass = 1.0;


	double ks = 7.0;
	double kd = 1.0;

    double structural_ks = ks;
    double structural_kd = kd;
    double flexion_ks = ks;
    double flexion_kd = kd;
	double shear_ks = ks;
    double shear_kd = kd;

	double deformation_rate = 1.1;

	ClothParticleSystem* system = new ClothParticleSystem(posX, posY, width, height, dist, mass ,structural_ks, structural_kd, flexion_ks, flexion_kd, shear_ks, shear_kd, deformation_rate);

	system->fixPoint(0,0);
	system->fixPoint(9,0);

	system->setIntegrationHook(RungeKutta);
	system->setDt(0.02);
	system->addForce(new GravityForce(system->getParticles(), 0.05));

	return system;
}

ParticleSystem* cloth2() {
	double posX = -0.5;
    double posY = 0.5;
    int width = 10;
    int height = 10;
    double dist = 0.1;
    double mass = 1.0;


	double ks = 7.0;
	double kd = 1.0;

    double structural_ks = ks;
    double structural_kd = kd;
    double flexion_ks = ks;
    double flexion_kd = kd;
	double shear_ks = ks;
    double shear_kd = kd;

	double deformation_rate = 1.1;

	ClothParticleSystem* system = new ClothParticleSystem(posX, posY, width, height, dist, mass ,structural_ks, structural_kd, flexion_ks, flexion_kd, shear_ks, shear_kd, deformation_rate);

	system->fixPointToHorizontalLine(0,0);
	system->fixPointToHorizontalLine(9,0);

	system->setIntegrationHook(RungeKutta);
	system->setDt(0.02);
	system->addForce(new GravityForce(system->getParticles(), 0.05));

	return system;
}

ParticleSystem* cloth3() {
	double posX = -0.5;
    double posY = 0.5;
    int width = 10;
    int height = 10;
    double dist = 0.1;
    double mass = 1.0;


	double ks = 7.0;
	double kd = 1.0;

    double structural_ks = ks;
    double structural_kd = kd;
    double flexion_ks = ks;
    double flexion_kd = kd;
	double shear_ks = ks;
    double shear_kd = kd;

	double deformation_rate = 1.1;

	ClothParticleSystem* system = new ClothParticleSystem(posX, posY, width, height, dist, mass ,structural_ks, structural_kd, flexion_ks, flexion_kd, shear_ks, shear_kd, deformation_rate);

	system->fixPoint(0,0);
	
	system->addWall(new Wall(Vec2f(-0.65, 0), 0.5*PI, 0.1));

	system->setIntegrationHook(RungeKutta);
	system->setDt(0.02);
	system->addForce(new GravityForce(system->getParticles(), 0.05));

	return system;
}

ParticleSystem* rigid1(float dt) {
	ParticleSystem* system = new ParticleSystem();

	system->setDt(dt);

	// Set integration scheme.
	system->setIntegrationHook(Euler);


	system->addRigid(new myRectangle(0.9,0.1,0.05,Vec2f(0.5,0.5)));

	// Create Global Constraint Jacobian Matrix
	Constraint::GlobalJ = new GlobalMatrix(0, system->particleCount() * 2);
	Constraint::GlobalJdot = new GlobalMatrix(0, system->particleCount() * 2);
	Constraint::global_cons_num = 0;
	Constraint::kd = 0.2;
	Constraint::ks = 0.3;

	return system;
}