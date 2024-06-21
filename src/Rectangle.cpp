#include "Rectangle.h"

myRectangle::myRectangle(double mass, float width, float height, Vec2f ConstructPos) : RigidBody(ConstructPos) {	
	Mass = mass;
	m_width = width;
	m_height = height;
	m_torque = 0;

	m_Vertices.push_back(new Particle(Vec2f(-m_width / 2, -m_height / 2), mass/4.0));
	m_Vertices.push_back(new Particle(Vec2f(m_width / 2, -m_height / 2), mass/4.0));
	m_Vertices.push_back(new Particle(Vec2f(m_width / 2, m_height / 2), mass/4.0));
	m_Vertices.push_back(new Particle(Vec2f(-m_width / 2, m_height / 2), mass/4.0));

	initInertia();
	reset();
}

void myRectangle::initInertia() {
	m_Inertia = Mass * (m_width * m_width + m_height * m_height) / 12;
}
