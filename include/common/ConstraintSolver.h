#pragma once

#include "Particle.h"
#include "Constraint.h"

#include <vector>

void ApplyConstraintForce(std::vector<Particle*> pVector, std::vector<Constraint*> cVector);
void GetGlobalVectors(std::vector<Particle*> pVector, double* q, double* Q, double* W, double* qdot);
void GetGlobalMatrices(std::vector<Constraint*> cVector, double* C, double* Cdot);
double calculateError(int n, double* x, double* b, double* W);