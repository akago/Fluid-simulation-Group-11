#include "Force.h"

double degreesToRadians(double degrees) {
	return degrees * (M_PI / 180.0);
}

AngularSpring::AngularSpring(Particle* x1, Particle* x2, Particle* x3, double angle, double ks, double kd) :
	endPoint1(x1), anglePoint(x2), endPoint2(x3), m_angle(angle), m_ks(ks), m_kd(kd) {}

Vec2f AngularSpring::getMidpoint() {
	return Vec2f((endPoint1->m_Position[0] + endPoint2->m_Position[0]) / 2, (endPoint1->m_Position[1] + endPoint2->m_Position[1]) / 2);
}

double AngularSpring::_distance(const Vec2f &x1, const Vec2f &x2) {
	return sqrt(pow(x1[0] - x2[0], 2) + pow(x1[1] - x2[1], 2));
}

double AngularSpring::calculate_rest_len() {
	double x2x1 = _distance(anglePoint->m_Position, endPoint1->m_Position);
	double x2x3 = _distance(anglePoint->m_Position, endPoint2->m_Position);

	//printf("cos 60: %.8f\n", cos(m_angle));
	return sqrt((x2x1 * x2x1) + (x2x3 * x2x3) + 2 * (x2x1 * x2x3) * cos(m_angle)) / 2;

}

Vec2f AngularSpring::getUnitVector(const Vec2f &x1, const Vec2f &x2) {
	double x1x2 = _distance(x1, x2);
	return Vec2f((x1[0] - x2[0]) / x1x2, (x1[1] - x2[1]) / x1x2);
}

double AngularSpring::applyForce() {
	double d0;  //Rest length between angle point and midpoint
	double F; // Spring force on endpoint
	Vec2f uv = getUnitVector(endPoint1->m_Position, endPoint2->m_Position); // Normalised vector in direction of endpoint1->endpoint2
	Vec2f midposition = getMidpoint();
	Particle* mid = new Particle(midposition);
	mid->reset(); // Set m_Position

	d0 = calculate_rest_len();
	SpringForce* virtual_force = new SpringForce(mid, anglePoint, d0, m_ks, m_kd);

	double d = _distance(mid->m_Position, anglePoint->m_Position);
	/*printf("midpoint: (%.2f,%.2f)\n", mid->m_Position[0], mid->m_Position[1]);
	printf("p1: (%.2f,%.2f)\n", endPoint1->m_Position[0], endPoint1->m_Position[1]);
	printf("p3: (%.2f,%.2f)\n", endPoint2->m_Position[0], endPoint2->m_Position[1]);
	printf("��AngularSpring��==================Current distance��%.8f, Current rest distance: %.8f\n============", d, d0);*/
	F = virtual_force->applyForce() / 2;

	endPoint1->m_Force[0] -= uv[0] * F;
	endPoint1->m_Force[1] -= uv[1] * F;
	endPoint2->m_Force[0] += uv[0] * F;
	endPoint2->m_Force[1] += uv[1] * F;

	free(mid);
	free(virtual_force);

	return 0.0;
}

void AngularSpring::draw() {

}
