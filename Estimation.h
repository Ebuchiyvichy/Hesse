#pragma once
#include "Condition.h"

//find octahedral norme of vector
float	octah_vect_norm(float *x, int size)
{
	float	norm = 0;

	for (int i = 0; i < size; i++)
		norm += fabs(x[i]);
	return (norm);
}

//find cube norme of vector
float	cube_vect_norm(float *x, int size)
{
	float	norm = 0;

	for (int i = 0; i < size; i++)
		if (norm < x[i])
			norm = x[i];
	return (norm);
}

//find sferical norm of vector
float	sfer_vect_norm(float *x, int DIM)
{
	float	norm = 0;

	for (int i = 0; i < DIM; i++)
		norm += x[i] * x[i];
	return (sqrt(norm));
}

//find absolutely difference of two vectors
float	*abs_diff_vector(float *a, float *b, int DIM)
{
	float	*c = new float[DIM];

	for (int i = 0; i < DIM; i++)
		c[i] = fabs(a[i] - b[i]);
	return (c);
}

//find difference of two vectors
float	*diff_vector(float *a, float *b, int DIM)
{
	float	*c = new float[DIM];

	for (int i = 0; i < DIM; i++)
		c[i] = a[i] - b[i];
	return (c);
}

//find conditional number's estimate

//for cube norme
float	cube_estimate_number(float *b, const Matrix& R, const Matrix& T)
{
	float	pert = 0.1;
	float	*b_, *tmp = new float[T.size];
	float	*x, *x_ = NULL;
	float	*dx = NULL;
	float	my_dx = 0;
	float	*db;
	float	nbr;

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
	nbr = (my_dx / cube_vect_norm(x, T.size)) / (cube_vect_norm(db, T.size) / cube_vect_norm(b, T.size));
	delete[] x;
	if (x_)
		delete[] x_;
	if (dx)
		delete[] dx;
	delete[] db;
	delete[] tmp;
	return (nbr);
}

//for octahedral norme
float	octah_estimate_number(float *b, const Matrix& R, const Matrix& T)
{
	float	pert = 0.1;
	int DIM = T.size;
	float	*b_, *tmp = new float[DIM];
	float	*x, *x_ = NULL;
	float	*dx = NULL;
	float	my_dx = 0;
	float	*db;
	float	nbr;

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
	nbr = (my_dx / octah_vect_norm(x, T.size)) / (octah_vect_norm(db, T.size) / octah_vect_norm(b, T.size));
	delete[] tmp;
	delete[] x;
	if (x_)
		delete[] x_;
	if (dx)
		delete[] dx;
	delete[] db;
	return (nbr);
}

//for sferical norm
float	sfer_estimate_number(float *b, const Matrix& R, const Matrix& T)
{
	float	pert = 0.1;
	int DIM = T.size;
	float	*b_, *tmp = new float[DIM];
	float	*x, *x_ = NULL;
	float	*dx = NULL;
	float	my_dx = 0;
	float	*db;
	float	nbr;

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
	nbr = (my_dx / sfer_vect_norm(x, T.size)) / (sfer_vect_norm(db, T.size) / sfer_vect_norm(b, T.size));
	delete[] tmp;
	delete[] x;
	if (x_)
		delete[] x_;
	if (dx)
		delete[] dx;
	delete[] db;
	return (nbr);
}