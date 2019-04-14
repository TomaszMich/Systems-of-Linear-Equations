#pragma once
#include "matrix.h"

class Equation
{
private:
	Matrix A, b, x;
	unsigned int N;

public:
	Equation(unsigned int N, int a1, int a2, int a3);
	void Jacobi();
	void GaussSeidel();
	void LUFactorization();
	void printX();
};