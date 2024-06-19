#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include <math.h> 
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <gfx/vec2.h>
#include<vector>

// Karen's CGD

#define MAX_STEPS 100


// Matrix class the solver will accept
class implicitMatrix
{
 public:
  virtual void matVecMult(double x[], double r[]) = 0;
};

// Matrix class the solver will accept
class implicitMatrixWithTrans : public implicitMatrix
{
 public:
  virtual void matVecMult(double x[], double r[]) = 0;
  virtual void matTransVecMult(double x[], double r[]) = 0;
};



class GlobalMatrix : public implicitMatrixWithTrans 
{
private:
	struct SparseBlock {
		int i;
		int j;
		int ilength;
		int jlength;
		std::vector<double> values;

		SparseBlock(int i_val, int j_val, std::vector<double> new_values) : i(i_val), j(j_val), ilength(1), jlength(2), values(new_values) {}
	};
	std::vector<SparseBlock> sparse_mat;
	int m_row;
	int m_col;
	
public:
	GlobalMatrix();
	GlobalMatrix(int row_num, int col_num);
	void updateRownum();
	void addEmptyBlock(int i, int j);
	void fillBlockat(int i, int j, std::vector<double> values);
	void matVecMult(double x[], double r[]) override;
	void matTransVecMult(double x[], double r[]) override;
	void printMatrix();
};


// Solve Ax = b for a symmetric, positive definite matrix A
// A is represented implicitely by the function "matVecMult"
// which performs a matrix vector multiple Ax and places result in r
// "n" is the length of the vectors x and b
// "epsilon" is the error tolerance
// "steps", as passed, is the maximum number of steps, or 0 (implying MAX_STEPS)
// Upon completion, "steps" contains the number of iterations taken
double ConjGrad(int n, implicitMatrix *A, double x[], double b[], 
		double epsilon,	// how low should we go?
		int    *steps);

// additional parameter: W - mass matrix
//						 m - length of W
double Constraint_ConjGrad(int n, int m, GlobalMatrix *A, double x[], double b[], double W[],
	double epsilon,	// how low should we go?
	int    steps);


// Some vector helper functions
float crossProduct(const Vec2f a, const Vec2f b);
void vecAddEqualWithFactor(int n, double r[], double v[], double factor);
void vecAddEqual(int n, double r[], double v[]);
void vecDiffEqual(int n, double r[], double v[]);
void vecAssign(int n, double v1[], double v2[]);
void vecTimesScalar(int n, double v[], double s);
void vecElewiseProd(int n, double r[], double v[]);
double vecDot(int n, double v1[], double v2[]);
double vecSqrLen(int n, double v[]);
void vecAddEqualWithFactor(int n, double r[], double v[], double factor);


#endif
