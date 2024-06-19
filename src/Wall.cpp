#include "Wall.h"

#define EPSILON 0.001

Wall::Wall(Vec2f position, double angle, double r) :
    position(position), normal(Vec2f(sin(angle), cos(angle))), r(r)
{
    double a = normal[0];
    double b = normal[1];
    double c = a*position[0]+b*position[1];
    
    int index = 0;

    if(b != 0) {
        double y = (c-a)/b;
        if(y >= -1 && y <= 1) {
            borderIntersection[index] = 1;
            borderIntersection[index+1] = y;
            index += 2;
        }
        y = (c+a)/b;
        if(y >= -1 && y <= 1 && borderIntersection) {
            borderIntersection[index] = -1;
            borderIntersection[index+1] = y;
            index += 2;
        }
    }

    if(a != 0) {
        double x = (c-b)/a;
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

void Wall::adjustParticle(Particle* particle, Vec2f prevPosition) {
    double c = (particle->m_Position-position)*normal;
    if(c <= 0) {
        if((prevPosition[0] != particle->m_Position[0] || prevPosition[1] != particle->m_Position[1]) && c < -EPSILON) {
            double lambda = (normal[0]*(position[0]-particle->m_Position[0]) + normal[1]*(position[1]-particle->m_Position[1])) / (normal[0]*(prevPosition[0]-particle->m_Position[0]) + normal[1]*(prevPosition[1]-particle->m_Position[1]));
            particle->m_Position = particle->m_Position + lambda * (prevPosition-particle->m_Position);
        }
    
        Vec2f x_n = (particle->m_Velocity*normal)*normal;
        Vec2f x_t = particle->m_Velocity - x_n;
        particle->m_Velocity = -r * x_n + x_t;
    }
}

void Wall::draw()
{
	glBegin( GL_LINES );
    glColor3f(0.0,0.0,1.0);
    glVertex2f( borderIntersection[0], borderIntersection[1] );
    glColor3f(0.0,0.0,1.0); 
    glVertex2f( borderIntersection[2], borderIntersection[3] );
    glEnd();
}


