#include "RigidBody.h"
#include "linearSolver.h"
#include <GL/glut.h>


RigidBody::RigidBody(Vec2f ConstructPos) 
	: m_ConstructPos(ConstructPos), m_Force(Vec2f(0.0,0.0)), m_Angle(0.0), m_AngularMomentum(0.0), m_LinearMomentum(Vec2f(0.0,0.0)), m_torque(0.0)
{
	updateVectors();
}

void RigidBody::updateVectors() {
	matrix R = matrix(Vec2f(cos(m_Angle), -sin(m_Angle)), Vec2f(sin(m_Angle), cos(m_Angle)));

	for (auto vertex : m_Vertices) {
		vertex->reset();
		vertex->m_Position = (R.multiByVec2(vertex->m_Position)) + m_Position;
		vertex->m_Velocity[0] = m_LinearMomentum[0] / Mass;
		vertex->m_Velocity[1] = m_LinearMomentum[1] / Mass;
	}
}

void RigidBody::clearForce() {
	m_Force = Vec2f(0.0, 0.0);
	m_torque = 0.0;
	updateVectors();
}

void RigidBody::draw() {
	glColor3f(0.0f, 1.0f, 0.0);
	glBegin(GL_POLYGON);
	for (auto vertex : m_Vertices) {
		glVertex2f(vertex->m_Position[0], vertex->m_Position[1]);
	}
	glEnd();
}


std::vector<float> RigidBody::getBBox() {
	std::vector<float> minmax = { 3.0,3.0,-3.0,-3.0 };
	
	for (auto vertex : m_Vertices) {
		if (minmax[0] > vertex->m_Position[0]) minmax[0] = vertex->m_Position[0]; 
		if (minmax[1] > vertex->m_Position[1]) minmax[1] = vertex->m_Position[1]; 
		if (minmax[2] < vertex->m_Position[0]) minmax[2] = vertex->m_Position[0];
		if (minmax[3] < vertex->m_Position[1]) minmax[3] = vertex->m_Position[1]; 
	}

	return minmax;
}


void RigidBody::addForce() {
	for (auto vertex : m_Vertices) {
		m_Force += vertex->m_Force;
	}
}

void RigidBody::addTorque() {
	for (auto vertex : m_Vertices) {
		m_torque += crossProduct(vertex->m_Position-m_Position, vertex->m_Force);
	}
}

void RigidBody::reset() {
	m_Force = Vec2f(0.0, 0.0);
	m_torque = 0.0;
	m_Angle = 0.0;
	m_AngularMomentum = 0.0;
	m_LinearMomentum = Vec2f(0.0, 0.0);

	m_Position = m_ConstructPos;

	updateVectors();
}


bool RigidBody::verticesInCell(const Vec2f& cell_center, const float h) {

	for (auto vertex : m_Vertices) {
		Vec2f point = vertex->m_Position;

		if (point[0] >= cell_center[0] - h * 0.5 && point[0] <= cell_center[0] + h * 0.5 &&
			point[1] >= cell_center[1] - h * 0.5 && point[1] <= cell_center[1] + h * 0.5)
		{
			return true;
		}
	}
	return false;
}


bool RigidBody::pointInPolygon(const Vec2f& point, const float h) {
	int intersections = 0;
	
	for (int i = 0; i < m_Vertices.size(); i++) {
		Vec2f v1 = m_Vertices[i]->m_Position;
		Vec2f v2 = m_Vertices[(i + 1) % m_Vertices.size()]->m_Position;

		if (((v1[1] <= point[1] && point[1] < v2[1]) || (v2[1] <= point[1] && point[1] < v1[1])) &&
			point[0] < (v2[0] - v1[0]) * (point[1] - v1[1]) / (v2[1] - v1[1]) + v1[0]) {
			intersections++;
		}
	}
	return (intersections % 2 != 0);
}

bool RigidBody::cellVerticesInPolygon(const Vec2f& cell_center, const float h) {
	
	return (pointInPolygon(Vec2f(cell_center[0] - h * 0.5, cell_center[1] - h * 0.5), h) ||
		pointInPolygon(Vec2f(cell_center[0] + h * 0.5, cell_center[1] - h * 0.5), h) ||
		pointInPolygon(Vec2f(cell_center[0] - h * 0.5, cell_center[1] + h * 0.5), h) ||
		pointInPolygon(Vec2f(cell_center[0] + h * 0.5, cell_center[1] + h * 0.5), h));
}

void RigidBody::projectToGrid(int N, int* internal_bd, float* bnd_vel_u, float* bnd_vel_v) {
	
	float h = 1.0 / N;
	std::vector<float> bbox = getBBox();

	/*printf("��projectToGrid�� h=%f\n",h);*/
	
	/* Alignment to grid cell offset */
	bbox[0] = floor(bbox[0] / h);
	bbox[1] = floor(bbox[1] / h);
	bbox[2] = ceil(bbox[2] / h);
	bbox[3] = ceil(bbox[3] / h);

	/* AABB Detection */
	/*printf("��projectToGrid�� start AABB Detection");
	printf("��projectToGrid�� From (%d,%d) to (%d,%d)\n", bbox[0], bbox[1], bbox[2], bbox[3]);*/
	// iterate grid cells overlapping with [minx,maxx] 
	for (int cell_x = bbox[0]; cell_x < bbox[2]; cell_x++) { 
		// iterate grid cells overlapping with [miny,maxy] 
		for (int cell_y = bbox[1]; cell_y < bbox[3]; cell_y++) { 
			Vec2f cell_center(cell_x * h + h * 0.5, cell_y * h + h * 0.5);
			if (cellVerticesInPolygon(cell_center, h) || verticesInCell(cell_center, h)) {
				// Vec2f vel = m_Velocity + (cell_center - m_Position) * m_Angular_Vec;
				bnd_vel_u[IX(cell_x, cell_y)] = m_LinearMomentum[0]/Mass;
				bnd_vel_v[IX(cell_x, cell_y)] = m_LinearMomentum[1]/Mass;
				internal_bd[IX(cell_x, cell_y)] = 2;
			}			
		}
	}
}



