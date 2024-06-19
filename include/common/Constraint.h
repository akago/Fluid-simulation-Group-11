#pragma once

#include "Particle.h"
#include "LinearSolver.h"
#include <cmath>

class Constraint {
public:
	static GlobalMatrix* GlobalJ;
	static GlobalMatrix* GlobalJdot;
	static int global_cons_num;
	static double ks;
	static double kd;

	virtual double eval() = 0;
	virtual double getTimeDeriv() = 0;
	virtual void fillJacobBlock() = 0;
	virtual void fillJacobDotBlock() = 0;
	virtual void draw() = 0;	 
};

/*
	C(x,y) = (x-xc)^2 + (y-yc)^2-r^2
	Cdot = 2(x-xc)vx + 2(y-yc)vy
*/
class CircularWireConstraint : public Constraint {
 public:
  CircularWireConstraint(int p_idx, Particle *p, const Vec2f & center, const double radius);
  double eval() override;
  double getTimeDeriv() override;
  void fillJacobBlock() override;
  void fillJacobDotBlock() override;
  void draw() override;

 private:
  int m_p_idx;
  Particle * const m_p;
  Vec2f const m_center;
  double const m_radius;
  // index of constraint instance 
  int m_c_idx;
};

/*
	Constraint Equation: C(x,y) = (a*x) + (b*y) - c
*/
class LineWireConstraint : public Constraint {
 public:
  LineWireConstraint(int p_idx, Particle *p, const double a, const double b, const double c);
  double eval() override;
  double getTimeDeriv() override;
  void fillJacobBlock() override;
  void fillJacobDotBlock() override;
  void draw() override;

 private:
  int m_p_idx;
  Particle * const m_p;
  const double m_a;
  const double m_b;
  const double m_c;
  // index of constraint instance 
  int m_c_idx;

  double borderIntersection[4];
};

/*
	Constraint Equation: C(x1,y1,x2,y2) = (x1-x2)^2 + (y1-y2)^2 - r^2
*/
class RodConstraint : public Constraint {
public:
	RodConstraint(int p1_idx, int p2_idx, Particle *p1, Particle * p2, double dist);
	double eval() override;
	double getTimeDeriv() override;
	void fillJacobBlock() override;
	void fillJacobDotBlock() override;
	void draw();

protected:
	int m_p1_idx;
	int m_p2_idx;
	Particle * const m_p1;
	Particle * const m_p2;
	double const m_dist;
	// index of constraint instance 
	int m_c_idx;
};

/*
	Constraint Equation: C(x1,y1,x2,y2) = ��[(x1-x2)^2 + (y1-y2)^2] - r
*/
class RodConstraintv2 : public RodConstraint {
public:
	using RodConstraint::RodConstraint;
	double eval() override;
	double getTimeDeriv() override;
	void fillJacobBlock() override;
	void fillJacobDotBlock() override;
};

