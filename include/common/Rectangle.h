#include "RigidBody.h"

class myRectangle : public RigidBody{
public:
	myRectangle(double mass, float width, float height, Vec2f ConstructPos);
	void initInertia();

protected:
	float m_width;
	float m_height;
};

