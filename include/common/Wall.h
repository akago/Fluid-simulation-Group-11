#pragma once

#include "Particle.h"
#include <gfx/vec2.h>
#include <GL/glut.h>
#include <cmath>

class Wall {
    public:
        Wall(Vec2f position, double angle, double r);
        void draw();
        void adjustParticle(Particle* particle, Vec2f prevPosition);
    private:
        const double r;
        const Vec2f position;
        const Vec2f normal;
        double borderIntersection[4];
};