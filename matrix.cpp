#include <math.h>
#include "matrix.h"

Matrix::Matrix()
{
	this->firstDim = 1;
	this->secondDim = 1;
	allocMemory();
}

Matrix::Matrix(const Matrix& m)
:Matrix(m.firstDim, m.secondDim)
{
	for (int i = 0; i < firstDim; i++)
	{
		for (int j = 0; j < secondDim; j++)
			M[i][j] = m[i][j];
	}
}

Matrix::Matrix(unsigned int firstDim, unsigned int secondDim)
{
	this->firstDim = firstDim;
	this->secondDim = secondDim;
	allocMemory();
}

void Matrix::allocMemory()
{
	this->M = new double*[firstDim];
	for (int i = 0; i < firstDim; ++i)
		this->M[i] = new double[secondDim];
}

void Matrix::createA(int a1, int a2, int a3)
{
	for (int i = 0; i < firstDim; i++)
	{
		for (int j = 0; j < secondDim; j++)
		{
			if (i == j)
				M[i][j] = a1;
			else if (abs(i - j) == 1)
				M[i][j] = a2;
			else if (abs(i - j) == 2)
				M[i][j] = a3;
			else
				M[i][j] = 0;
		}
	}
}

void Matrix::createB()
{
	for (int n = 0; n < firstDim; n++)
		M[n][0] = sin(n * 2);
}

void Matrix::createX()
{
	for (int n = 0; n < firstDim; n++)
		M[n][0] = 1;
}

double Matrix::norm(Matrix m)
{
	double value = 0;
	for (int i = 0; i < m.getFirstDim(); i++)
	{
		value += pow(m[i][0], 2);
	}
	return sqrt(value);
}

double* Matrix::operator[](size_t el)
{
	return M[el];
}

double* Matrix::operator[](size_t el)const
{
	return M[el];
}

Matrix Matrix::operator*(Matrix &m)
{
	Matrix temp(m);
	for (int i = 0; i < this->firstDim; i++)
		for (int j = 0; j < m.getSecondDim(); j++)
			temp[i][j] = 0;

	for (int i = 0; i < this->firstDim; ++i)
		for (int j = 0; j < m.getSecondDim(); ++j)
			for (int k = 0; k < this->secondDim; ++k)
				temp[i][j] += this->M[i][k] * m[k][j];

	return temp;
}

Matrix Matrix::operator-(Matrix &m)
{
	for (int i = 0; i < firstDim; i++)
		for (int j = 0; j < secondDim; j++)
			this->M[i][j] -= m[i][j];

	return *this;
}

Matrix Matrix::operator=(Matrix m)
{
	if (this == &m)
		return *this;
	this->~Matrix();
	this->firstDim = m.getFirstDim();
	this->secondDim = m.getSecondDim();
	allocMemory();
	for (int i = 0; i < firstDim; i++)
	{
		for (int j = 0; j < secondDim; j++)
			M[i][j] = m[i][j];
	}
	return *this;
}

unsigned int Matrix::getFirstDim()
{
	return firstDim;
}

unsigned int Matrix::getSecondDim()
{
	return secondDim;
}

Matrix::~Matrix()
{
	for (int i = 0; i < firstDim; i++)
		delete[] M[i];

	delete[] M;
}