#pragma once

class Matrix
{
private:
	unsigned int firstDim, secondDim;
	double** M;

public:
	Matrix(unsigned int firstDim, unsigned int secondDim);
	Matrix(const Matrix& m);
	Matrix();
	~Matrix();
	void allocMemory();
	void createA(int a1, int a2, int a3);
	void createB();
	void createX();
	static double norm(Matrix res);
	double* operator[](size_t el);
	double* operator[](size_t el) const;
	Matrix operator*(Matrix& m);
	Matrix operator-(Matrix& m);
	Matrix operator=(Matrix m);
	unsigned int getFirstDim();
	unsigned int getSecondDim();
};