#include <stdio.h>
#include <chrono>
#include "equation.h"

#define e 1e-9
#define LIMIT 1000

Equation::Equation(unsigned int N, int a1, int a2, int a3)
{
	this->N = N;
	A = Matrix(N, N);
	b = Matrix(N, 1);
	x = Matrix(N, 1);

	A.createA(a1, a2, a3);
	b.createB();
	x.createX();
}

void Equation::Jacobi()
{
	int itCount = 0;
	double value;
	Matrix xNew(N, 1);
	xNew.createX();
	auto start = std::chrono::high_resolution_clock::now();
	while (itCount < LIMIT)
	{
		for (int i = 0; i < N; i++)
		{
			value = b[i][0];
			for (int j = 0; j < N; j++)
			{
				if (j != i)
					value -= A[i][j] * x[j][0];
			}
			value /= A[i][i];
			xNew[i][0] = value;
		}
		x = xNew;
		itCount++;
		
		if (Matrix::norm((A * x) - b) <= e)
		{
			std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
			printf("Jacobi method took %d iterations and %fs.\n", itCount, elapsed.count());
			break;
		}
	}
}

void Equation::GaussSeidel()
{
	int itCount = 0;
	double value;
	auto start = std::chrono::high_resolution_clock::now();
	while (itCount < LIMIT)
	{
		for (int i = 0; i < N; i++)
		{
			value = b[i][0];
			for (int j = 0; j < N; j++)
			{
				if (j != i)
					value -= A[i][j] * x[j][0];
			}
			value /= A[i][i];
			x[i][0] = value;
		}
		itCount++;

		if (Matrix::norm((A * x) - b) <= e)
		{
			std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
			printf("Gauss-Seider method took %d iterations and %fs.\n", itCount, elapsed.count());
			break;
		}
	}
}

void Equation::LUFactorization()
{
	// decomposition of matrix
	Matrix lu(N, N);
	double sum = 0;
	for (int i = 0; i < N; i++)
	{
		for (int j = i; j < N; j++)
		{
			sum = 0;
			for (int k = 0; k < i; k++)
				sum += lu[i][k] * lu[k][j];
			lu[i][j] = A[i][j] - sum;
		}
		for (int j = i + 1; j < N; j++)
		{
			sum = 0;
			for (int k = 0; k < i; k++)
				sum += lu[j][k] * lu[k][i];
			lu[j][i] = (1 / lu[i][i]) * (A[j][i] - sum);
		}
	}

	// lu = L+U-I
	// find solution of Ly = b
	Matrix y(N, 1);
	for (int i = 0; i < N; i++)
	{
		sum = 0;
		for (int k = 0; k < i; k++)
			sum += lu[i][k] * y[k][0];
		y[i][0] = b[i][0] - sum;
	}
	// find solution of Ux = y
	for (int i = N - 1; i >= 0; i--)
	{
		sum = 0;
		for (int k = i + 1; k < N; k++)
			sum += lu[i][k] * x[k][0];
		x[i][0] = (1 / lu[i][i]) * (y[i][0] - sum);
	}
}

void Equation::printX()
{
	printf("x = [");
	for (int i = 0; i < N - 1; i++)
		printf("%f, ", x[i][0]);
	printf("%f]\n", x[N - 1][0]);
	printf("Norm: %E\n", Matrix::norm((A * x) - b));
}