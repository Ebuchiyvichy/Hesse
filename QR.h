#pragma once

#include "MatrixClass.h"

double	*cpy_vector(double *tmp, double *x, int size)
{

	for (int i = 0; i < size; i++)
		tmp[i] = x[i];
	return (tmp);
}

bool	degenerate_matrix(const Matrix& R)
{
	for (int i = 0; i < R.size; i++)
		if (fabs(R.value[R.mymap[i] - R.value + i]) < EPS)
			return (true);
	return (false);
}

// reverse Hause
double	*find_x(const Matrix& A, double *b_)
{
	double	*x = new double[A.size];

	for (int i = A.size - 1; i >= 0; i--)
	{
		for (int j = i + 1; j < A.size; j++)
			b_[i] = b_[i] - A.value[A.mymap[i] - A.value + j] * b_[j];
		b_[i] = b_[i] / A.value[A.mymap[i] - A.value + i];
	}
	x = cpy_vector(x, b_, A.size);
	return (x);
}

// find matrix T
void Matrix::QR_find_x(Matrix& A)
{
	double	c;
	double	s;
	double	tmp;

	for (int k = 0; k < size; k++)
	{
		for (int i = k + 1; i < size; i++)
		{
			if (fabs(A.value[A.mymap[k] - A.value + i]) > 10e-8)
			{
				c = A.value[A.mymap[k] - A.value + k] / sqrt(A.value[A.mymap[k] - A.value + k] * A.value[A.mymap[k] - A.value + k] + A.value[A.mymap[i] - A.value + k] * A.value[A.mymap[i] - A.value + k]);
				s = A.value[A.mymap[i] - A.value + k] / sqrt(A.value[A.mymap[k] - A.value + k] * A.value[A.mymap[k] - A.value + k] + A.value[A.mymap[i] - A.value + k] * A.value[A.mymap[i] - A.value + k]);
				for (int j = 0; j <= i; j++) // change T-matrix
				{
					tmp = value[mymap[k] - value + j];
					value[mymap[k] - value + j] = value[mymap[k] - value + j] * c + value[mymap[i] - value + j] * s;
					value[mymap[i] - value + j] = c * value[mymap[i] - value + j] - s * tmp;
				}
				for (int j = k; j < size; j++) // change A-matrix
				{
					tmp = A.value[A.mymap[k] - A.value + j];
					A.value[A.mymap[k] - A.value + j] = c * A.value[A.mymap[k] - A.value + j] + s * A.value[A.mymap[i] - A.value + j];
					A.value[A.mymap[i] - A.value + j] = c * A.value[A.mymap[i] - A.value + j] - s * tmp;
				}
			}
		}
	}
	std::cout << "QR method matrix R:" << std::endl;
	A.print();

	std::cout << "QR method matrix T:" << std::endl;
	print();
}