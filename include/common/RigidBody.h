#pragma once

#include <gfx/vec2.h>
#include "Matrix.h"
#include "Particle.h"
#include <vector>

#define IX(i,j) ((i)+(N+2)*(j))


class RigidBody
{
public:
	RigidBody(Vec2f ConstructPos);

	void draw();
	void clearForce();
	std::vector<float> getBBox();
	void updateR(float angle);
	void projectToGrid(int N, int* internal_bd, float* bnd_vel_u, float* bnd_vel_v);
	bool verticesInCell(const Vec2f& cell_center, const float h);
	bool pointInPolygon(const Vec2f& point, const float h);
	bool cellVerticesInPolygon(const Vec2f& cell_center, const float h);
	void addForce();
	void addTorque();
	void reset();

	/* Vertices of the rigid*/
	std::vector<Particle*> m_Vertices;

	/* Constant quantities */
	double Mass;
	Vec2f m_ConstructPos; /* x0 */
	float m_Inertia; /* I */

	/* State variables */
	Vec2f m_Position;	 /* x(t) */
	matrix m_R;			/* R(t) */	
	Vec2f m_Velocity; /* v(t) */
	float m_Angular_Vec; /* ω(t) */

	/* Computed quantities */
	Vec2f m_Force; //(right, up), same as opengl coordinate
	float m_torque;
	float m_Angle; // the angle that a rigidbody has been rotated since last timestep

};


