#include "Force.h"
#include <GL/glut.h>

MouseForce::MouseForce(Vec2f mousePos, std::vector<Particle*> particles, double dist, double ks, double kd) :
	m_mousePos(mousePos), m_pVector(particles), m_dist(dist), m_ks(ks), m_kd(kd)
{
	leftMouseDown = false;
	rightMouseDown = false;
}


void MouseForce::draw()
{

}

double MouseForce::applyForce()
{
	if (leftMouseDown)
	{
		for (int i = 0; i < m_pVector.size(); i++)
		{
			Vec2f dist = m_pVector[i]->m_Position - m_mousePos;
			Vec2f vecDif = m_pVector[i]->m_Velocity;

			float distance = sqrt(dist[0] * dist[0] + dist[1] * dist[1]);
			float dotProduct = dist * vecDif;

			float scalar = (m_ks * (distance - m_dist) + m_kd * (dotProduct / distance));
			Vec2f result = scalar * (dist / distance);

			m_pVector[i]->m_Force -= result;
		}

	}
	else if (rightMouseDown)
	{
		for (int i = 0; i < m_rigid_bodies.size(); i++)
		{
			Vec2f dist = m_rigid_bodies[i]->m_Position - m_mousePos;
			Vec2f vecDif = m_rigid_bodies[i]->m_Velocity;
																				
			float distance = sqrt(dist[0] * dist[0] + dist[1] * dist[1]);
			float dotProduct = dist * vecDif;

			float scalar = (m_ks * (distance - m_dist) + m_kd * (dotProduct / distance));

			Vec2f result = scalar * (dist / distance);

			m_rigid_bodies[i]->m_Force -= result;
		}
	}
	return 0;
}

void MouseForce::newMousePosition(Vec2f mousePos)
{
	m_mousePos = mousePos;
}

void MouseForce::selectParticles(std::vector<Particle*> pVector)
{
	m_pVector = pVector;
	leftMouseDown = true;
}

void MouseForce::clearParticles()
{
	leftMouseDown = false;
}

void MouseForce::selectRigidbodies(std::vector<RigidBody*> rigid_bodies)
{
	m_rigid_bodies = rigid_bodies;
	rightMouseDown = true;
}

void MouseForce::clearRigidbodies()
{
	rightMouseDown = false;
}