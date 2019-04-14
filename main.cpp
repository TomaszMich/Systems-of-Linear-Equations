#include "equation.h"

//Implementation of iterative (Jacobi and Gauss-Seidel) and direct LU factorization (decomposition) methods of solving systems of linear equations.

int main()
{
	const unsigned int N = 990;
	int a1 = 13, a2 = -1, a3 = -1;
	
	Equation equation(N, a1, a2, a3);
	equation.Jacobi();
	//equation.GaussSeidel();
	//equation.LUFactorization();

	equation.printX();

	return 0;
}

