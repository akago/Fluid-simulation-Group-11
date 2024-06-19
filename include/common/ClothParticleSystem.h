#pragma once

#include "ParticleSystem.h"
#include "Particle.h"
#include "Force.h"
#include "Constraint.h"

#include <vector>
#include <gfx/vec2.h>

const double TWOSQUAREROOT = 1.414213562373095;

class ClothParticleSystem : public ParticleSystem {
public:
    ClothParticleSystem(double posX, double posY, int width, int height, double dist, double mass, double structural_ks, double structural_kd, double flexion_ks, double flexion_kd, double shear_ks, double shear_kd, double deformation_rate);
    ~ClothParticleSystem() override;
    void fixPoint(int x, int y);
    void fixPointToHorizontalLine(int x, int y);
    void fixPointToVerticalLine(int x, int y);
    void fixPointToLine(int x, int y, double angle);

    void simulationStep() override;
private:
    double const posX;
    double const posY;
    int const width;
    int const height;
    double const dist;
    double const deformationRate;

    std::vector<bool> isConstrained;

    int index(int x, int y);
    Vec2f pos(int x, int y);

    void AdjustElongatedParticles(int p_idx0, int p_idx1, double maxDistance);
};