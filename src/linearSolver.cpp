#include "linearSolver.h"
#include "Constraint.h"
#include "util.h"

// vector helper functions

float crossProduct(const Vec2f a, const Vec2f b)
{
	return a[0] * b[1] - a[1] * b[0];
}

void vecAddEqualWithFactor(int n, double r[], double v[], double factor)
{
  for (int i = 0; i < n; i++)
    r[i] = r[i] + factor * v[i];
}

void vecAddEqual(int n, double r[], double v[])
{
  for (int i = 0; i < n; i++)
    r[i] = r[i] + v[i];
}

void vecDiffEqual(int n, double r[], double v[])
{
  for (int i = 0; i < n; i++)
    r[i] = r[i] - v[i];
}

void vecAssign(int n, double v1[], double v2[])
{
  for (int i = 0; i < n; i++)
    v1[i] = v2[i];
}

void vecTimesScalar(int n, double v[], double s)
{
  for (int i = 0; i < n; i++)
    v[i] *= s;
}

void vecElewiseProd(int n, double r[], double v[]) {
	for (int i = 0; i < n; i++)
		r[i] = r[i] * v[i];
}

double vecDot(int n, double v1[], double v2[])
{
  double dot = 0;
  for (int i = 0; i < n; i++)
    dot += v1[i] * v2[i];
  return dot;
}

double vecSqrLen(int n, double v[])
{
  return vecDot(n, v, v);
}

double ConjGrad(int n, implicitMatrix *A, double x[], double b[], 
		double epsilon,	// how low should we go?
		int    *steps)
{
  int		i, iMax;
  double	alpha, beta, rSqrLen, rSqrLenOld, u;

  double *r = (double *) malloc(sizeof(double) * n);
  double *d = (double *) malloc(sizeof(double) * n);
  double *t = (double *) malloc(sizeof(double) * n);
  double *temp = (double *) malloc(sizeof(double) * n);

  vecAssign(n, x, b);

  vecAssign(n, r, b);

  A->matVecMult(x, temp);
  
  vecDiffEqual(n, r, temp);

  rSqrLen = vecSqrLen(n, r);

  vecAssign(n, d, r);

  i = 0;
  if (*steps)
    iMax = *steps;
  else
    iMax = MAX_STEPS;
		
  if (rSqrLen > epsilon)
    while (i < iMax) {	
      i++;
      A->matVecMult(d, t);
      u = vecDot(n, d, t);
      
      if (u == 0) {
	//printf("(SolveConjGrad) d'Ad = 0\n");
	break;
      }
      
      // How far should we go?
      alpha = rSqrLen / u;
      
      // Take a step along direction d
      vecAssign(n, temp, d);
      vecTimesScalar(n, temp, alpha);
      vecAddEqual(n, x, temp);
      
      if (i & 0x3F) {
	vecAssign(n, temp, t);
	vecTimesScalar(n, temp, alpha);
	vecDiffEqual(n, r, temp);
      } else {
	// For stability, correct r every 64th iteration
	vecAssign(n, r, b);
	A->matVecMult(x, temp);
	vecDiffEqual(n, r, temp);
      }
      
      rSqrLenOld = rSqrLen;
      rSqrLen = vecSqrLen(n, r);
      
      // Converged! Let's get out of here
      if (rSqrLen <= epsilon)
	break;			    
      
      // Change direction: d = r + beta * d
      beta = rSqrLen/rSqrLenOld;
      vecTimesScalar(n, d, beta);
      vecAddEqual(n, d, r);
    }
  
  // free memory

  free(r);
  free(d);
  free(t);
  free(temp);
		
  *steps = i;
  return(rSqrLen);
}



double Constraint_ConjGrad(int n, int m, GlobalMatrix *A, double x[], double b[], double W[],
	double epsilon,	// how low should we go?
	int    steps)
{
	int		i, iMax;
	double	alpha, beta, rSqrLen, rSqrLenOld, u;

	double *r = (double *)malloc(sizeof(double) * n);
	double *d = (double *)malloc(sizeof(double) * n);
	double *t = (double *)malloc(sizeof(double) * n);
	double *temp = (double *)malloc(sizeof(double) * n);
	// dimension of W: m
	double *temp1 = (double *)malloc(sizeof(double) * m);

	vecAssign(n, x, b);
	vecAssign(n, r, b);
	////printf("��Constraint_ConjGrad�� b = ");
	//printVector(n, b);
	A->matTransVecMult(x, temp1); 
	////printf("��Constraint_ConjGrad�� Jt x = ");
	//printVector(m, temp1);
	vecElewiseProd(m, temp1, W);
	////printf("��Constraint_ConjGrad�� W Jt x = ");
	//printVector(m, temp1);
	A->matVecMult(temp1, temp);
	////printf("��Constraint_ConjGrad�� J W Jt x = ");
	//printVector(n, temp);

	vecDiffEqual(n, r, temp);
	////printf("��Constraint_ConjGrad�� r = ");
	//printVector(n, r);
	rSqrLen = vecSqrLen(n, r);
	////printf("��Constraint_ConjGrad��rsqlen: %.8f, epsilon: %.8f\n", rSqrLen, epsilon);
	vecAssign(n, d, r);

	i = 0;
	if (steps)
		iMax = steps;
	else
		iMax = MAX_STEPS;
	if (rSqrLen > epsilon)
		while (i < iMax) {
			i++;
			////printf("��Constraint_ConjGrad��%d-th steps, error is %.6f\n", i, rSqrLen);
			//A->matVecMult(d, t);
			A->matTransVecMult(d, temp1);
			vecElewiseProd(m, temp1, W);
			A->matVecMult(temp1, t);

			u = vecDot(n, d, t);

			if (u == 0) {
				//printf("(SolveConjGrad) d'Ad = 0\n");
				break;
			}

			// How far should we go?
			alpha = rSqrLen / u;

			// Take a step along direction d
			vecAssign(n, temp, d);
			vecTimesScalar(n, temp, alpha);
			vecAddEqual(n, x, temp);

			if (i & 0x3F) {
				vecAssign(n, temp, t);
				vecTimesScalar(n, temp, alpha);
				vecDiffEqual(n, r, temp);
			}
			else {
				// For stability, correct r every 64th iteration
				vecAssign(n, r, b);
				//A->matVecMult(x, temp);
				A->matTransVecMult(x, temp1);
				vecElewiseProd(m, temp1, W);
				A->matVecMult(temp1, temp);
				vecDiffEqual(n, r, temp);
			}

			rSqrLenOld = rSqrLen;
			rSqrLen = vecSqrLen(n, r);

			// Converged! Let's get out of here
			if (rSqrLen <= epsilon)
				break;

			// Change direction: d = r + beta * d
			beta = rSqrLen / rSqrLenOld;
			vecTimesScalar(n, d, beta);
			vecAddEqual(n, d, r);
		}
	// free memory

	free(r);
	free(d);
	free(t);
	free(temp);
	free(temp1);
	//*steps = i;

	return(rSqrLen);
}

GlobalMatrix::GlobalMatrix() :
	m_row(0), m_col(0){}

GlobalMatrix::GlobalMatrix(int row_num, int col_num) :
	m_row(row_num), m_col(col_num) {}

void GlobalMatrix::addEmptyBlock(int i, int j) {
	sparse_mat.push_back(SparseBlock(i, j, std::vector<double>(2, 0.0)));
	updateRownum();
}

void GlobalMatrix::fillBlockat(int i, int j, std::vector<double> values) {
	for (int ii = 0; ii < sparse_mat.size(); ii++) {
		if (sparse_mat[ii].i == i && sparse_mat[ii].j == j) {
			sparse_mat[ii].values = values;
		}
	}
}

void GlobalMatrix::updateRownum() {
	m_row = Constraint::global_cons_num;
}

void GlobalMatrix::matVecMult(double x[], double r[]) {
	for (int ii = 0; ii < m_row; ii++) {
		r[ii] = 0.0;
	}
	for (auto block : sparse_mat) {
		for (int ii = 0; ii < block.jlength; ii++) {
			r[block.i] += block.values[ii] * x[block.jlength * block.j + ii];
		}
	}
}

void GlobalMatrix::matTransVecMult(double x[], double r[]) {
	for (int ii = 0; ii < m_col; ii++) {
		r[ii] = 0.0;
	}
	for (auto block : sparse_mat) {
		r[block.jlength * block.j] += block.values[0] * x[block.i];
		r[block.jlength * block.j + 1] +=  block.values[1] * x[block.i];
	}
}

void GlobalMatrix::printMatrix() {
	for (int ii = 0; ii < m_row; ii++) {
		for (int jj = 0; jj < m_col/2; jj++) {
			int flag = 0;
			for (auto block : sparse_mat) {
				if (block.i == ii && block.j == jj) {
					flag = 1;
					printf("%.8f  ", block.values[0]);
					printf("%.8f  ", block.values[1]);
				}
			}
			if (flag) {
				flag = 0;
			}
			else {
				double zero = 0.0;
				printf("%.8f  ", zero);
				printf("%.8f  ", zero);
			}
		}
		printf("\n");
	}
	printf("\n");
}