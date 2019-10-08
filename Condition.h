#include "QR.h"

//multiply matrix * vector
double	*multi_vect(double* I, const Matrix& T)
{
	double *tmp;

	tmp = new double[T.size];
	for (int i = 0; i < T.size; i++)
		tmp[i] = 0;
	for (int j = 0; j < T.size; j++)
		for (int i = 0; i < T.size; i++)
			tmp[j] += T.value[T.mymap[j]-T.value+i] * I[i];
	return (tmp);
}

//find inverse matrix
 void Matrix::inverse_matrix( const Matrix& R, Matrix& T)
{
	double *x;
	double *y;
	x = new double [T.size];
	y = new double [T.size];

	T.Tranc();
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
				std::cout << "Vector x in QR:" << std::endl;
				for (int j = 0; j < T.size; j++)
				{
					std::cout << std::setw(8) << x[j] << std::endl;
				}
				std::cout << std::endl;
		for (int j = 0; j < T.size; j++)
		{
			value[mymap[i] - value + j] = x[j];
		}
	}
	T.Tranc();
}

 //find cube norm
double	cube_norm(const Matrix& A)
{
	double	norm = 0;
	double	sum = 0;

	for (int i = 0; i < A.size; i++)
	{
		for (int j = 0; j < A.size; j++)
			sum += fabs(A.value[A.mymap[i]-A.value +j]);
		if (sum > norm)
			norm = sum;
		sum = 0;
	}
	return (norm);
}
 
//find octahedral norm
double	octah_norm(const Matrix& A)
{
	double	norm = 0;
	double	sum = 0;

	for (int i = 0; i < A.size; i++)
	{
		for (int j = 0; j < A.size; j++)
			sum += fabs(A.value[A.mymap[i]-A.value +j]);
		if (sum > norm)
			norm = sum;
		sum = 0;
	}
	return (norm);
}

//find sferical norm
double	sfer_norm(const Matrix& A)
{
	double	sum = 0;

	for (int i = 0; i < A.size; i++)
	{
		for (int j = 0; j < A.size; j++)
			sum += fabs(A.value[A.mymap[i]-A.value +j]);
	}
	return (sqrt(sum));
}
