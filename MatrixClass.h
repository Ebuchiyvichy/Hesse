#pragma once
//#define float float
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
float	EPS = 10e-8;

class Matrix
{
public:
	float	*value;
	float	**mymap;
	float	*rvalue;
	int		size;

public:
	Matrix()
	{
		std::ifstream file;
		file.open("matrix.txt");
		file >> size;
		value = new float[size*size];
		mymap = new float*[size];
		rvalue = new float[size];
		float number;
		int n = 0;

		while (n != size)
		{
			for (int i = 0; i < size; i++)
			{
				file >> number;
				value[i + n * size] = number;
			}
			file >> number;
			rvalue[n] = number;
			n++;
		}
		for (int i = 0; i != size; i++)
		{
			mymap[i] = &(value[size*i]);
		}
		file.close();
	}

	Matrix(const Matrix& A)
	{
		size = A.size;
		value = new float[size*size];
		mymap = new float*[size];
		rvalue = new float[size];
		for (int i = 0; i != (size*size); i++)
		{
			value[i] = A.value[i];
		}
		for (int j = 0; j != size; j++)
		{
			mymap[j] = &(value[j*size]);
			rvalue[j] = rvalue[j];
		}
	}

	Matrix(int s, int i = 1)
	{
		size = s;
		value = new float[size*size];
		mymap = new float*[size];
		rvalue = NULL;
		for (int i = 0; i != (size*size); i++)
		{
			if ((i % size) == (i / size))
				value[i] = 1;
			else
				value[i] = 0;
		}
		for (int i = 0; i != size; i++)
		{
			mymap[i] = &(value[i*size]);
		}
	}

	~Matrix()
	{
		delete[] value;
		delete[] mymap;
		delete[] rvalue;
	}

	void			Tranc();
	void			print();

	friend bool		degenerate_matrix(const Matrix& R);
	friend float	*multi_vect(float* I, const Matrix& T);

	//functions for QR
	void			QR_find_x(Matrix& A);
	void			inverse_matrix(const Matrix& R, Matrix& T);
	friend float	*find_x(const Matrix& A, float *b_);

	//functions for norm
	friend float	cube_norm(const Matrix& A);
	friend float	octah_norm(const Matrix& A);
	friend float	sfer_norm(const Matrix& A);

	//functions for cond nbr
	friend float	octah_estimate_number(float *b, const Matrix& R, const Matrix& T);
	friend float	cube_estimate_number(float *b, const Matrix& R, const Matrix& T);
	friend float	sfer_estimate_number(float *b, const Matrix& R, const Matrix& T);

	//functions for Gauss
	friend void		GaussLeft(Matrix& A);
	friend void		delim(Matrix &A, int k);
	friend void		vych(Matrix &A, int k);
	friend void		GaussRight(Matrix &A);
};

float	*cpy_vector(float *tmp, float *x, int size);
float	*abs_diff_vector(float *a, float *b, int DIM);

//functions for vector norm
float	octah_vect_norm(float *x, int size);
float	cube_vect_norm(float *x, int size);
float	sfer_vect_norm(float *x, int DIM);

void Matrix::print()
{
	for (int n = 0; n != size; n++)
	{
		for (int j = 0; j != size; j++)
		{
			if (fabs(value[mymap[n] - value + j]) < EPS)
				std::cout << std::setw(8) << 0 << "\t";
			else
				std::cout << std::setw(8) << value[mymap[n] - value + j] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void Matrix::Tranc()
{ // i - string, j - column
	float temp;

	for (int i = 0; i != size; i++)
	{
		for (int j = i + 1; j != size; j++)
		{
			temp = value[mymap[i] - value + j];
			value[(mymap[i] - value) + j] = value[(mymap[j] - value) + i];
			value[(mymap[j] - value) + i] = temp;
		}
	}
}

void delim(Matrix &A, int k)
{
	int i = k;
	float temp = A.value[A.mymap[i] - A.value + k];
	for (int j = k; j != A.size; j++)
	{
		A.value[A.mymap[i] - A.value + j] = A.value[A.mymap[i] - A.value + j] / temp;
	}
	A.rvalue[i] = A.rvalue[i] / temp;
}

void vych(Matrix &A, int k)
{
	for (int i = (k + 1); i != A.size; i++)
	{
		//columns
		float temp = A.value[A.mymap[i] - A.value + k];
		for (int j = 0; j != A.size; j++)
			A.value[A.mymap[i] - A.value + j] = A.value[A.mymap[i] - A.value + j] - temp * A.value[A.mymap[k] - A.value + j];
		A.rvalue[i] -= temp * A.rvalue[k];
	}
}

void GaussRight(Matrix &A)
{
	// for k-string our matrix
	for (int k = 0; k != A.size; k++)
	{
		//search max in column
		float max = A.value[A.mymap[k] - A.value + k];
		int maxstring = k;
		for (int i = k; i < A.size; i++)
		{
			if (A.value[A.mymap[i] - A.value + k] > max)
			{
				max = A.value[A.mymap[i] - A.value + k];
				maxstring = i;
			}
		}
		//change rvalue
		float temp = A.rvalue[k];
		A.rvalue[k] = A.rvalue[maxstring];
		A.rvalue[maxstring] = temp;
		// change strings
		int temparr = (A.mymap[k] - A.value);
		A.mymap[k] = &(A.value[A.mymap[maxstring] - A.value]);
		A.mymap[maxstring] = &(A.value[temparr]);
		delim(A, k);
		vych(A, k);
	}
}

void GaussLeft(Matrix &A)
{
	for (int i = A.size - 1; i >= 0; i--)
	{
		for (int j = i + 1; j < A.size; j++)
			A.rvalue[i] = A.rvalue[i] - A.value[A.mymap[i] - A.value + j] * A.rvalue[j];
		A.rvalue[i] = A.rvalue[i] / A.value[A.mymap[i] - A.value + i];
	}
}

void reversematrixgauss(Matrix &A, Matrix &E)
{
	for (int i = 0; i < A.size; i++)
	{
		Matrix A_cpy(A);
		for (int j = 0; j < A.size; j++)
		{
			if (j == i)
				A_cpy.rvalue[j] = 1;
			else
				A_cpy.rvalue[j] = 0;
		}
			GaussRight(A_cpy);
			GaussLeft(A_cpy);
		for (int j = 0; j < A.size; j++)
		{
			E.value[E.mymap[i] - E.value + j] = A_cpy.rvalue[j];
		}
	}
	E.Tranc();
}
