#pragma once
#include "QR.h"
#include "MatrixClass.h"

//multiply matrix * vector
float	*multi_vect(float* I, const Matrix& T)
{
	float *tmp;

	tmp = new float[T.size];
	for (int i = 0; i < T.size; i++)
		tmp[i] = 0;
	for (int j = 0; j < T.size; j++)
		for (int i = 0; i < T.size; i++)
			tmp[j] += T.value[T.mymap[j] - T.value + i] * I[i];
	return (tmp);
}

//find inverse matrix
void Matrix::inverse_matrix(const Matrix& R, Matrix& T)
{
	//with Hausse
	float *x;
	float *y;
	x = new float [T.size];
	y = new float [T.size];
	Matrix E(R.size, 1);

	for (int i = 0; i < T.size; i++)
	{
		for (int j = 0; j < T.size; j++)
		{
			if (j == i)
				y[j] = 1;
			else
				y[j] = 0;
		}
		x = find_x(R, multi_vect(y, T));
		for (int j = 0; j < T.size; j++)
		{
			value[mymap[i] - value + j] = x[j];
		}
	}
	Tranc();
	delete[] x;
	delete[] y;
/*	//with inverse triangle matrix
	Matrix E(R.size, 1);
	for (int k = R.size - 1; k >= 0; k--)
	{
		E.value[E.mymap[k] - E.value + k] = E.value[E.mymap[k] - E.value + k] / R.value[R.mymap[k] - R.value + k];
		for (int j = k - 1; j >= 0; j--)
			{
				for (int i = R.size - 1; i > j; i--)
					E.value[E.mymap[j] - E.value + k] -= R.value[R.mymap[j] - R.value + i] * E.value[E.mymap[i] - E.value + k];
				E.value[E.mymap[j] - E.value + k] /= R.value[R.mymap[j] - R.value + j];
			}	
	}
	//for (int k = 0; k < R.size; k++)
	//{
	//	value[mymap[k] - value + k] = value[mymap[k] - value + k] / R.value[R.mymap[k] - R.value + k];
	//	for (int j = k - 1; j >= 0; j--)
	//	{
	//		for (int i = R.size - 1; i > j; i--)
	//			value[mymap[j] - value + k] -= R.value[R.mymap[j] - R.value + i] * value[mymap[i] - value + k];
	//		value[mymap[j] - value + k] /= R.value[R.mymap[j] - R.value + j];
	//	}
	//}
	for (int i = 0; i < R.size; i++)
	{
		for (int j = 0; j < R.size; j++)
		{
			value[mymap[i] - value + j] = 0;
			for (int k = 0; k < R.size; k++)
				value[mymap[i] - value + j] = E.value[E.mymap[i] - E.value + k] * T.value[T.mymap[k] - T.value + j];
		}
	}*/
}

//find cube norm
float	cube_norm(const Matrix& A)
{
	float	norm = 0;
	float	sum = 0;

	for (int i = 0; i < A.size; i++)
	{
		for (int j = 0; j < A.size; j++)
			sum += fabs(A.value[A.mymap[i] - A.value + j]);
		if (sum > norm)
			norm = sum;
		sum = 0;
	}
	return (norm);
}

//find octahedral norm
float	octah_norm(const Matrix& A)
{
	float	norm = 0;
	float	sum = 0;

	for (int i = 0; i < A.size; i++)
	{
		for (int j = 0; j < A.size; j++)
			sum += fabs(A.value[A.mymap[j] - A.value + i]);
		if (sum > norm)
			norm = sum;
		sum = 0;
	}
	return (norm);
}

//find sferical norm
float	sfer_norm(const Matrix& A)
{
	float	sum = 0;

	for (int i = 0; i < A.size; i++)
	{
		for (int j = 0; j < A.size; j++)
			sum += fabs(A.value[A.mymap[i] - A.value + j]);
	}
	return (sqrt(sum));
}
