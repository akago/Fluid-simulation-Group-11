#include "ConstraintSolver.h"
#include "util.h"

void GetGlobalVectors(std::vector<Particle*> pVector, double* q, double* Q, double* W, double* qdot) {
	for (auto particle : pVector) {
		*(q++) = particle->m_Position[0];
		*(q++) = particle->m_Position[1];
		*(Q++) = particle->m_Force[0];
		*(Q++) = particle->m_Force[1];
		*(W++) = particle->mass;
		*(W++) = particle->mass;
		*(qdot++) = particle->m_Velocity[0];
		*(qdot++) = particle->m_Velocity[1];
	}
}

void GetGlobalMatrices(std::vector<Constraint*> cVector, double* C, double* Cdot) {
	for (int ii = 0; ii < cVector.size(); ii++) {
		cVector[ii]->fillJacobBlock();
		cVector[ii]->fillJacobDotBlock();
		*(C++) = cVector[ii]->eval();
		*(Cdot++) = cVector[ii]->getTimeDeriv();
	}
}

double calculateError(int n, double* x, double* b, double* W) {
	double* temp = (double*)malloc(n * sizeof(double));
	double* res = (double*)malloc(Constraint::global_cons_num * sizeof(double));
	Constraint::GlobalJ->matTransVecMult(x, temp);
	vecElewiseProd(n, temp, W);
	Constraint::GlobalJ->matVecMult(temp, res);

	printf("��calculateError��==========START===============\n");
	printf("JWJtx is: ");
	printVector(Constraint::global_cons_num, res);
	printf("b is: ");
	printVector(Constraint::global_cons_num, b);
	printf("��calculateError��==========END===============\n");

	free(temp);
	free(res);
}

void ApplyConstraintForce(std::vector<Particle*> pVector, std::vector<Constraint*> cVector) {
	double epsilon = 0.000000001;
	int steps = 0;
	int dimension = pVector.size() * 2;

	double* q = (double*)malloc(dimension * sizeof(double)); // global position vector
	double* Q = (double*)malloc(dimension * sizeof(double)); // global force vector
	double* W = (double*)malloc(dimension * sizeof(double)); // mass diagnoal matrix
	double* qdot = (double*)malloc(dimension * sizeof(double)); // global velocity vector

	double* C = (double*)malloc(Constraint::global_cons_num * sizeof(double));
	double* Cdot = (double*)malloc(Constraint::global_cons_num * sizeof(double));
	double* b = (double*)malloc(Constraint::global_cons_num * sizeof(double));
	double* x = (double*)malloc(Constraint::global_cons_num * sizeof(double));

	GetGlobalVectors(pVector, q, Q, W, qdot);

	// GetGlobalMatrix
	GetGlobalMatrices(cVector, C, Cdot);

	Constraint::GlobalJdot->matVecMult(qdot, b);						// Jdot*qdot
	vecElewiseProd(dimension, Q, W);									// W*Q
	Constraint::GlobalJ->matVecMult(Q, x);								// JWQ
	vecTimesScalar(Constraint::global_cons_num, C, Constraint::ks);		// ks*C
	vecTimesScalar(Constraint::global_cons_num, Cdot, Constraint::kd);  // kd*Cdot

	// Sum up
	vecAddEqual(Constraint::global_cons_num, b, x);
	vecAddEqual(Constraint::global_cons_num, b, C);
	vecAddEqual(Constraint::global_cons_num, b, Cdot);
	vecTimesScalar(Constraint::global_cons_num, b, -1.0);

	Constraint_ConjGrad(Constraint::global_cons_num, dimension, Constraint::GlobalJ, x, b, W, epsilon, steps);
	//calculateError(dimension, x, b, W);

	Constraint::GlobalJ->matTransVecMult(x, Q);
	double* constraint_forces = Q;
	for (auto particle : pVector) {
		particle->m_Force[0] += *(constraint_forces++);
		particle->m_Force[1] += *(constraint_forces++);
	}
	free(q);
	free(Q);
	free(W);
	free(qdot);
	free(C);
	free(Cdot);
	free(b);
	free(x);
}