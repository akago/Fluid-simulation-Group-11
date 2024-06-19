#include "Force.h"
#include "linearSolver.h"
#include <cmath>
#include <GL/glut.h>

SpringForce::SpringForce(Particle *p1, Particle *p2, double dist, double ks, double kd) :
  m_p1(p1), m_p2(p2), m_dist(dist), m_ks(ks), m_kd(kd) {}



double SpringForce::applyForce() {
	double *I = (double *)malloc(sizeof(double) * 2);
	double *Idot = (double *)malloc(sizeof(double) * 2);
	double *force = (double *)malloc(sizeof(double) * 2);
	double distance = 0.0;
		
	I[0] = m_p1->m_Position[0] - m_p2->m_Position[0];
	I[1] = m_p1->m_Position[1] - m_p2->m_Position[1];

	distance = sqrt(I[0] * I[0] + I[1] * I[1]);

	Idot[0] = m_p1->m_Velocity[0] - m_p2->m_Velocity[0];
	Idot[1] = m_p1->m_Velocity[1] - m_p2->m_Velocity[1];

	// Idot = ks(|I|-r) + Idot * I * kd / |I|
	force[0] = -(Idot[0] * I[0] * m_kd / distance + m_ks * (distance - m_dist)) * I[0] / distance;
	force[1] = -(Idot[1] * I[1] * m_kd / distance + m_ks * (distance - m_dist)) * I[1] / distance;

	// Apply force
	m_p1->m_Force[0] += force[0];
	m_p1->m_Force[1] += force[1];
	m_p2->m_Force[0] -= force[0];
	m_p2->m_Force[1] -= force[1];

	double forceMagnitude = sqrt((force[0] * force[0]) + (force[1] * force[1]));

	free(I);
	free(Idot);
	free(force);

	// Return magnitude of spring force
	return forceMagnitude;
}

void SpringForce::draw()
{
  glBegin( GL_LINES );
  glColor3f(0.6, 0.7, 0.8);
  glVertex2f( m_p1->m_Position[0], m_p1->m_Position[1] );
  glColor3f(0.6, 0.7, 0.8);
  glVertex2f( m_p2->m_Position[0], m_p2->m_Position[1] );
  glEnd();
}
