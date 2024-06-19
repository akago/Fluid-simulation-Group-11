#include "Force.h"

GravityForce::GravityForce(std::vector<Particle*> pVector, double gConstant) :
	m_pVector(pVector), m_gconstant(gConstant) {}

double GravityForce::applyForce() {
	for (auto p : m_pVector) {
		p->m_Force[1] -= p->mass * m_gconstant;
	}
	return m_gconstant;
}

void GravityForce::draw() {
	
}

