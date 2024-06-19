#include "Constraint.h"
#include <GL/glut.h>


RodConstraint::RodConstraint(int p1_idx, int p2_idx, Particle *p1, Particle * p2, double dist) :
 m_p1_idx(p1_idx), m_p2_idx(p2_idx), m_p1(p1), m_p2(p2), m_dist(dist), m_c_idx(global_cons_num++) 
{
	GlobalJ->addEmptyBlock(m_c_idx, m_p1_idx);
	GlobalJ->addEmptyBlock(m_c_idx, m_p2_idx);
	GlobalJdot->addEmptyBlock(m_c_idx, m_p1_idx);
	GlobalJdot->addEmptyBlock(m_c_idx, m_p2_idx);
}


double RodConstraint::eval() {
	// C(x1,y1,x2,y2) = (x1-x2)^2 + (y1-y2)^2 - r^2
	return pow(m_p1->m_Position[0] - m_p2->m_Position[0], 2) + pow(m_p1->m_Position[1] - m_p2->m_Position[1], 2) - pow(m_dist, 2);
}

double RodConstraint::getTimeDeriv() {
	// C' = dC/dt = 2*(x1-x2)*(vx1-vx2) + 2*(y1-y2)*(vy1-vy2)
	return 2 * (m_p1->m_Position[0] - m_p2->m_Position[0]) * (m_p1->m_Velocity[0] - m_p2->m_Velocity[0]) + 2 * (m_p1->m_Position[1] - m_p2->m_Position[1]) * (m_p1->m_Velocity[1] - m_p2->m_Velocity[1]);
}

void RodConstraint::fillJacobBlock() {
	// dC/d(q1) = (dC/d(x1), dC/d(y1)) = (2(x1-x2) ,2(y1-y2))
	// dC/d(q2) = (dC/d(x2), dC/d(y2)) = (-2(x1-x2) ,-2(y1-y2))
	std::vector<double> newval;

	newval.push_back((double)2 * (m_p1->m_Position[0] - m_p2->m_Position[0]));
	newval.push_back((double)2 * (m_p1->m_Position[1] - m_p2->m_Position[1]));
	GlobalJ->fillBlockat(m_c_idx, m_p1_idx, newval);

	newval[0] = -newval[0];
	newval[1] = -newval[1];
	GlobalJ->fillBlockat(m_c_idx, m_p2_idx, newval);
	//printf("¡¾fillJacobBlock¡¿ GlobalJ:\n");
	//GlobalJ->printMatrix();
}

void RodConstraint::fillJacobDotBlock() {
	// J' = dC'/d(q1) = (2(vx1-vx2), 2(vy1-vy2))
	// J' = dC'/d(q2) =  (-2(vx1-vx2), -2(vy1-vy2))
	std::vector<double> newval;

	newval.push_back((double)2 * (m_p1->m_Velocity[0] - m_p2->m_Velocity[0]));
	newval.push_back((double)2 * (m_p1->m_Velocity[1] - m_p2->m_Velocity[1]));
	GlobalJdot->fillBlockat(m_c_idx, m_p1_idx, newval);

	newval[0] = -newval[0];
	newval[1] = -newval[1];
	
	GlobalJdot->fillBlockat(m_c_idx, m_p2_idx, newval);
	//printf("¡¾fillJacobDotBlock¡¿ GlobalJdot:\n");
	//GlobalJdot->printMatrix();
}



void RodConstraint::draw()
{
  glBegin( GL_LINES );
  glColor3f(0.8, 0.7, 0.6);
  glVertex2f( m_p1->m_Position[0], m_p1->m_Position[1] );
  glColor3f(0.8, 0.7, 0.6);
  glVertex2f( m_p2->m_Position[0], m_p2->m_Position[1] );
  glEnd();

}


/*
	RodConstraintv2
*/


double RodConstraintv2::eval() {
	// C(x1,y1,x2,y2) = ¡Ì[(x1-x2)^2 + (y1-y2)^2] - r
	return sqrt(pow(m_p1->m_Position[0] - m_p2->m_Position[0], 2) + pow(m_p1->m_Position[1] - m_p2->m_Position[1], 2)) - m_dist;
}

double RodConstraintv2::getTimeDeriv() {
	// C' = dC/dt = ((x1-x2)(vx1-vx2) + (y1-y2)(vy1-vy2)) / (¡Ì[(x1-x2)^2 + (y1-y2)^2]) 
	double term1 = sqrt(pow(m_p1->m_Position[0] - m_p2->m_Position[0], 2) + pow(m_p1->m_Position[1] - m_p2->m_Position[1], 2));
	double term2 = (m_p1->m_Position[0] - m_p2->m_Position[0]) * (m_p1->m_Velocity[0] - m_p2->m_Velocity[0]) + (m_p1->m_Position[1] - m_p2->m_Position[1]) * (m_p1->m_Velocity[1] - m_p2->m_Velocity[1]);
	return term2 / term1;
}

void RodConstraintv2::fillJacobBlock() {
	// dC/d(q1) = (dC/d(x1), dC/d(y1)) = (x1-x2) / (¡Ì[(x1-x2)^2 + (y1-y2)^2])  , (y1-y2) / (¡Ì[(x1-x2)^2 + (y1-y2)^2)
	// dC/d(q2) = (dC/d(x2), dC/d(y2)) = -(x1-x2) / (¡Ì[(x1-x2)^2 + (y1-y2)^2])  , -(y1-y2) / (¡Ì[(x1-x2)^2 + (y1-y2)^2)
	std::vector<double> newval;

	double term1 = sqrt(pow(m_p1->m_Position[0] - m_p2->m_Position[0], 2) + pow(m_p1->m_Position[1] - m_p2->m_Position[1], 2));

	newval.push_back((double)2 * (m_p1->m_Position[0] - m_p2->m_Position[0]) / term1);
	newval.push_back((double)2 * (m_p1->m_Position[1] - m_p2->m_Position[1]) / term1);
	GlobalJ->fillBlockat(m_c_idx, m_p1_idx, newval);

	newval[0] = -newval[0];
	newval[1] = -newval[1];
	GlobalJ->fillBlockat(m_c_idx, m_p2_idx, newval);
}

void RodConstraintv2::fillJacobDotBlock() {
	// u = (x1-x2)^2 + (y1-y2)^2, u' = 2(x1-x2)(vx1-vx2) + 2(y1-y2)(vy1-vy2)
	// term = 1/¡Ìu, term' = (-1/2) * u^(-3/2) * u' = -((x1-x2)(vx1-vx2) + (y1-y2)(vy1-vy2)) / u^(3/2), 
	// J =  (x1-x2) * term1, J' = dC'/d(q1) = (term1' * (x1-x2) + term1 * (vx1-vx2), term1' * (y1-y2) + term1 * (vy1-vy2))
	// J' = dC'/d(q2) =  -(term1' * (x1-x2) + term1 * (vx1-vx2)), -(term1' * (y1-y2) + term1 * (vy1-vy2))
	std::vector<double> newval;

	//double u = pow(m_p1->m_Position[0] - m_p2->m_Position[0], 2) + pow(m_p1->m_Position[1] - m_p2->m_Position[1], 2);
	//double udot = 2 * (m_p1->m_Position[0] - m_p2->m_Position[0]) * (m_p1->m_Velocity[0] - m_p2->m_Velocity[0]) + 2 * (m_p1->m_Position[1] - m_p2->m_Position[1]) * (m_p1->m_Velocity[1] - m_p2->m_Velocity[1]);
	//double term1 = 1 / sqrt(u);
	//double term1deriv = -udot / (2 * pow(u, 3.0 / 2.0));
	//newval.push_back(term1deriv * (m_p1->m_Position[0] - m_p2->m_Position[0]) + term1 * m_p1->m_Velocity[0] - m_p2->m_Velocity[0]);
	//newval.push_back(term1deriv * (m_p1->m_Position[1] - m_p2->m_Position[1]) + term1 * m_p1->m_Velocity[1] - m_p2->m_Velocity[1]);
	double term1 = (m_p1->m_Velocity[0] - m_p2->m_Velocity[0]) * pow(m_p1->m_Position[1] - m_p2->m_Position[1], 2);
	double term2 = (m_p1->m_Position[0] - m_p2->m_Position[0]) * (m_p1->m_Position[1] - m_p2->m_Position[1]) * (m_p1->m_Velocity[1] - m_p2->m_Velocity[1]);
	double term3 = (m_p1->m_Velocity[1] - m_p2->m_Velocity[1]) * pow(m_p1->m_Position[0] - m_p2->m_Position[0], 2);
	double term4 = (m_p1->m_Position[1] - m_p2->m_Position[1]) * (m_p1->m_Position[0] - m_p2->m_Position[0]) * (m_p1->m_Velocity[0] - m_p2->m_Velocity[0]);
	double term5 = pow(pow(m_p1->m_Position[0] - m_p2->m_Position[0], 2) + pow(m_p1->m_Position[1] - m_p2->m_Position[1], 2), 1.5);
	newval.push_back((term1 - term2) / term5);
	newval.push_back((term3 - term4) / term5);
	
	GlobalJdot->fillBlockat(m_c_idx, m_p1_idx, newval);

	newval[0] = -newval[0];
	newval[1] = -newval[1];
	GlobalJdot->fillBlockat(m_c_idx, m_p2_idx, newval);
}
