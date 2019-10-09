#pragma once
//#define double float
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
double	EPS = 10e-8;

class Matrix
{
public:
	double	*value;
	double	**mymap;
	double	*rvalue;
	int		size;

public:
	Matrix()
	{
		std::ifstream file;
		file.open("matrix.txt");
		file >> size;
		value = new double[size*size];
		mymap = new double*[size];
		rvalue = new double[size];
		double number;
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
		value = new double[size*size];
		mymap = new double*[size];
		rvalue = new double[size];
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
		value = new double[size*size];
		mymap = new double*[size];
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
	friend double	*multi_vect(double* I, const Matrix& T);

	//functions for QR
	void			QR_find_x(Matrix& A);
	void			inverse_matrix(const Matrix& R, Matrix& T);
	friend double	*find_x(const Matrix& A, double *b_);

	//functions for norm
	friend double	cube_norm(const Matrix& A);
	friend double	octah_norm(const Matrix& A);
	friend double	sfer_norm(const Matrix& A);

	//functions for cond nbr
	friend double	octah_estimate_number(double *b, const Matrix& R, const Matrix& T);
	friend double	cube_estimate_number(double *b, const Matrix& R, const Matrix& T);
	friend double	sfer_estimate_number(double *b, const Matrix& R, const Matrix& T);

	//functions for Gauss
	friend void		GaussLeft(Matrix& A);
	friend void		delim(Matrix &A, int k);
	friend void		vych(Matrix &A, int k);
	friend void		GaussRight(Matrix &A);
};

double	*cpy_vector(double *tmp, double *x, int size);
double	*abs_diff_vector(double *a, double *b, int DIM);

//functions for vector norm
double	octah_vect_norm(double *x, int size);
double	cube_vect_norm(double *x, int size);
double	sfer_vect_norm(double *x, int DIM);

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
	double temp;

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
	double temp = A.value[A.mymap[i] - A.value + k];
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
		double temp = A.value[A.mymap[i] - A.value + k];
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
		double max = A.value[A.mymap[k] - A.value + k];
		int maxstring = k;
		for (int i = k; i != A.size; i++)
		{
			if (A.value[A.mymap[i] - A.value + k] > max)
			{
				max = A.value[A.mymap[i] - A.value + k];
				maxstring = i;
			}
		}
		//change rvalue
		double temp = A.rvalue[k];
		A.rvalue[k] = A.rvalue[maxstring];
		A.rvalue[maxstring] = temp;
		// change strings
		int temparr = (A.mymap[k] - A.value);
		A.mymap[k] = &(A.value[A.mymap[maxstring] - A.value]);
		A.mymap[maxstring] = &(A.value[temparr]);
		delim(A, k);
		vych(A, k);
	}
	std::cout << "Matrix A after straight run of Gauss:" << std::endl;
	A.print();
}

void GaussLeft(Matrix &A)
{
	for (int i = A.size - 1; i >= 0; i--)
	{
		for (int j = i + 1; j < A.size; j++)
			A.rvalue[i] = A.rvalue[i] - A.value[A.mymap[i] - A.value + j] * A.rvalue[j];
		A.rvalue[i] = A.rvalue[i] / A.value[A.mymap[i] - A.value + i];
	}
	std::cout << "Vector in Gauss x:" << std::endl;
	for (int i = 0; i != A.size; i++)
		std::cout << std::setw(8) << A.rvalue[i] << std::endl;
	std::cout << std::endl;
}
//doesn't work
void reversematrixgauss(Matrix &A, Matrix &E)
{
	for (int k = 0; k < A.size; k++)
	{
		//search max in column
		double max = A.value[A.mymap[k] - A.value + k];
		int maxstring = k;
		for (int i = k; i < A.size; i++)
		{
			if (A.value[A.mymap[i] - A.value + k] > max)
			{
				max = A.value[A.mymap[i] - A.value + k];
				maxstring = i;
			}
		}
		//change strings for matrix A
		int temparr = (A.mymap[k] - A.value);
		A.mymap[k] = &(A.value[A.mymap[maxstring] - A.value]);
		A.mymap[maxstring] = &(A.value[temparr]);
		//change for matrix E right 
		temparr = (E.mymap[k] - E.value);
		E.mymap[k] = &(E.value[E.mymap[maxstring] - E.value]);
		E.mymap[maxstring] = &(E.value[temparr]);
		//delim 
		int i = k;
		double temp = A.value[A.mymap[i] - A.value + k];
		for (int j = k; j < A.size; j++)
		{
			A.value[A.mymap[i] - A.value + j] = A.value[A.mymap[i] - A.value + j] / temp;
			E.value[E.mymap[i] - E.value + j] = E.value[E.mymap[i] - E.value + j] / temp;
		}
		//difference
		for (int i = (k + 1); i < A.size; i++)
		{
			//columns
			temp = A.value[A.mymap[i] - A.value + k];
			for (int j = 0; j < A.size; j++)
			{
				A.value[A.mymap[i] - A.value + j] -= temp * A.value[A.mymap[k] - A.value + j];
				E.value[E.mymap[i] - E.value + j] -= temp * E.value[E.mymap[k] - E.value + j];
			}
		}
	}
	
/*	for (int k = A.size - 1; k > 0; k--)
	{
		double	tmp1 = A.value[A.mymap[k] - A.value + k];
		E.value[E.mymap[k] - E.value + k] = E.value[E.mymap[k] - E.value + k] / tmp1;
		A.value[A.mymap[k] - A.value + k] = A.value[A.mymap[k] - A.value + k] / tmp1;
		for (int j = 0; j < k; j++)
		{
			double	tmp = A.value[A.mymap[j] - A.value + k] / tmp1;
			A.value[A.mymap[j] - A.value + k] = A.value[A.mymap[j] - A.value + k] - tmp * A.value[A.mymap[k] - A.value + k];
			E.value[E.mymap[j] - E.value + k] = E.value[E.mymap[j] - E.value + k] - tmp * E.value[E.mymap[k] - E.value + k];
		}
	}*/

	for (int k = A.size - 1; k >= 0; k--)
	{
		for (int j = k - 1; j >= 0; j--)
		{
			for (int i = A.size - 1; i > j; i--)
				E.value[E.mymap[j] - E.value + k] -= A.value[A.mymap[j] - A.value + i] * E.value[E.mymap[i] - E.value + k];
		}
	}
}
