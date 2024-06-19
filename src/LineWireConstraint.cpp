#include "Constraint.h"
#include <GL/glut.h>

LineWireConstraint::LineWireConstraint(int p_idx, Particle *p, const double a, const double b, const double c) :
	m_p_idx(p_idx), m_p(p),  m_a(a), m_b(b), m_c(c), m_c_idx(global_cons_num++)
{
	GlobalJ->addEmptyBlock(m_c_idx, m_p_idx);
	GlobalJdot->addEmptyBlock(m_c_idx, m_p_idx);

    // Compute intersections with window border
    int index = 0;

    if(m_b != 0) {
        double y = (c-m_a)/m_b;
        if(y >= -1 && y <= 1) {
            borderIntersection[index] = 1;
            borderIntersection[index+1] = y;
            index += 2;
        }
        y = (c+m_a)/m_b;
        if(y >= -1 && y <= 1 && borderIntersection) {
            borderIntersection[index] = -1;
            borderIntersection[index+1] = y;
            index += 2;
        }
    }

    if(m_a != 0) {
        double x = (c-m_b)/m_a;
        if(x > -1 && x < 1) {
            borderIntersection[index] = x;
            borderIntersection[index+1] = 1;
            index += 2;
        }
        if(x > -1 && x < 1) {
            borderIntersection[index] = x;
            borderIntersection[index+1] = -1;
            index += 2;
        }
    }
}

double LineWireConstraint::eval() {
	// C(x,y) = a*x + b*y - c
	return m_a*m_p->m_Position[0] + m_b*m_p->m_Position[1] - m_c;
}

double LineWireConstraint::getTimeDeriv() {
	// C' = dC/dt = a*vx + b*vy 
	return m_a*m_p->m_Velocity[0] + m_b*m_p->m_Velocity[1];
}

void LineWireConstraint::fillJacobBlock() {
	// compute (dC/dx, dC/dy) = (a , b)
	std::vector<double> newval;
	newval.push_back(m_a);
	newval.push_back(m_b);

	GlobalJ->fillBlockat(m_c_idx, m_p_idx, newval);
}

void LineWireConstraint::fillJacobDotBlock() {
	// J' = dC'/dq = (0, 0)
	std::vector<double> newval;
	newval.push_back(0);
	newval.push_back(0);

	GlobalJdot->fillBlockat(m_c_idx, m_p_idx, newval);
}


void LineWireConstraint::draw()
{
	glBegin( GL_LINES );
    glColor3f(0.0,1.0,0.0);
    glVertex2f( borderIntersection[0], borderIntersection[1] );
    glColor3f(0.0,1.0,0.0); 
    glVertex2f( borderIntersection[2], borderIntersection[3] );
    glEnd();
}


