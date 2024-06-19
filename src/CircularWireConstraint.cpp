#include "Constraint.h"
#include <GL/glut.h>

#define PI 3.1415926535897932384626433832795

static void draw_circle(const Vec2f & vect, float radius)
{
	glBegin(GL_LINE_LOOP);
	glColor3f(0.0,1.0,0.0); 
	for (int i=0; i<360; i=i+18)
	{
		float degInRad = i*PI/180;
		glVertex2f(vect[0]+cos(degInRad)*radius,vect[1]+sin(degInRad)*radius);
	}
	glEnd();
}

CircularWireConstraint::CircularWireConstraint(int p_idx, Particle *p, const Vec2f & center, const double radius) :
	m_p_idx(p_idx), m_p(p),  m_center(center), m_radius(radius), m_c_idx(global_cons_num++)
{
	GlobalJ->addEmptyBlock(m_c_idx, m_p_idx);
	GlobalJdot->addEmptyBlock(m_c_idx, m_p_idx);
}

double CircularWireConstraint::eval() {
	// C(x,y) = (x-xc)^2 + (y-yc)^2 - r^2
	return pow(m_p->m_Position[0] - m_center[0], 2) + pow(m_p->m_Position[1] - m_center[1], 2) - pow(m_radius, 2);
}

double CircularWireConstraint::getTimeDeriv() {
	// C' = dC/dt = 2(x-xc)*vx + 2(y-yc)*vy
	return 2 * (m_p->m_Position[0] - m_center[0]) * m_p->m_Velocity[0] + 2 * (m_p->m_Position[1] - m_center[1]) * m_p->m_Velocity[1];
}

void CircularWireConstraint::fillJacobBlock() {
	// compute (dC/dx, dC/dy) = (2(x-xc) ,2(y-yc))
	std::vector<double> newval;
	newval.push_back((double)2 * (m_p->m_Position[0] - m_center[0]));
	newval.push_back((double)2 * (m_p->m_Position[1] - m_center[1]));

	GlobalJ->fillBlockat(m_c_idx, m_p_idx, newval);
}

void CircularWireConstraint::fillJacobDotBlock() {
	// J' = dC'/dq = (2vx, 2vy)
	std::vector<double> newval;
	newval.push_back((double)2 * m_p->m_Velocity[0]);
	newval.push_back((double)2 * m_p->m_Velocity[1]);

	GlobalJdot->fillBlockat(m_c_idx, m_p_idx, newval);
}


void CircularWireConstraint::draw()
{
	draw_circle(m_center, m_radius);
}


