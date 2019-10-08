#pragma once
#include "Condition.h"

//find octahedral norme of vector
double	octah_vect_norm(double *x, int size)
{
	double	norm = 0;

	for (int i = 0; i < size; i++)
		norm += fabs(x[i]);
	return (norm);
}

//find cube norme of vector
double	cube_vect_norm(double *x, int size)
{
	double	norm = 0;

	for (int i = 0; i < size; i++)
		if (norm < x[i])
			norm = x[i];
	return (norm);
}

//find sferical norm of vector
double	sfer_vect_norm(double *x, int DIM)
{
	double	norm = 0;

	for (int i = 0; i < DIM; i++)
		norm += x[i] * x[i];
	return (sqrt(norm));
}

//find absolutely difference of two vectors
double	*abs_diff_vector(double *a, double *b, int DIM)
{
	double	*c = new double[DIM];

	for (int i = 0; i < DIM; i++)
		c[i] = fabs(a[i] - b[i]);
	return (c);
}

//find difference of two vectors
double	*diff_vector(double *a, double *b, int DIM)
{
	double	*c = new double[DIM];

	for (int i = 0; i < DIM; i++)
		c[i] = a[i] - b[i];
	return (c);
}

//find conditional number's estimate

//for cube norme
double	cube_estimate_number(double *b, const Matrix& R, const Matrix& T)
{
	double	pert = 0.1;
	double	*b_, *tmp = new double[T.size];
	double	*x, *x_;
	double	*dx;
	double	my_dx = 0;
	double	*db;

	b_ = multi_vect(b, T);
	x = find_x(R, b_);
	for (int i = 0; i < 2 * T.size; i++)
	{
		tmp = cpy_vector(tmp, b, T.size);
		if (i < T.size)
			tmp[i] += pert;
		else
			tmp[i - T.size] -= pert;
		b_ = multi_vect(tmp, T);
		x_ = find_x(R, b_);
		dx = abs_diff_vector(x, x_, T.size);
		if (my_dx < cube_vect_norm(dx, T.size))
			my_dx = cube_vect_norm(dx, T.size);
		delete[] b_;
	}
	db = abs_diff_vector(tmp, b, T.size);
	return ((my_dx / cube_vect_norm(x, T.size)) / (cube_vect_norm(db, T.size) / cube_vect_norm(b, T.size)));
}

//for octahedral norme
double	octah_estimate_number(double *b, const Matrix& R, const Matrix& T)
{
	double	pert = 0.1;
	int DIM = T.size;
	double	*b_, *tmp = new double[DIM];
	double	*x, *x_;
	double	*dx;
	double	my_dx = 0;
	double	*db;

	b_ = multi_vect(b, T);
	x = find_x(R, b_);
	for (int i = 0; i < 2 * DIM; i++)
	{
		tmp = cpy_vector(tmp, b, DIM);
		if (i < DIM)
			tmp[i] += pert;
		else
			tmp[i - DIM] -= pert;
		b_ = multi_vect(tmp, T);
		x_ = find_x(R, b_);
		dx = abs_diff_vector(x, x_, DIM);
		if (my_dx < octah_vect_norm(dx, DIM))
			my_dx = octah_vect_norm(dx, DIM);
	}
	db = abs_diff_vector(tmp, b, DIM);
	return ((my_dx / octah_vect_norm(x, DIM)) / (octah_vect_norm(db, DIM) / octah_vect_norm(b, DIM)));
}

//for sferical norm
double	sfer_estimate_number(double *b, const Matrix& R, const Matrix& T)
{
	double	pert = 0.1;
	int DIM = T.size;
	double	*b_, *tmp = new double[DIM];
	double	*x, *x_;
	double	*dx;
	double	my_dx = 0;
	double	*db;

	b_ = multi_vect(b, T);
	x = find_x(R, b_);
	for (int i = 0; i < 2 * DIM; i++)
	{
		tmp = cpy_vector(tmp, b, DIM);
		if (i < DIM)
			tmp[i] += pert;
		else
			tmp[i - DIM] -= pert;
		b_ = multi_vect(tmp, T);
		x_ = find_x(R, b_);
		dx = abs_diff_vector(x, x_, DIM);
		if (my_dx < sfer_vect_norm(dx, DIM))
			my_dx = sfer_vect_norm(dx, DIM);
		delete[] b_;
	}
	db = abs_diff_vector(tmp, b, DIM);
	return ((my_dx / sfer_vect_norm(x, DIM)) / (sfer_vect_norm(db, DIM) / sfer_vect_norm(b, DIM)));
}