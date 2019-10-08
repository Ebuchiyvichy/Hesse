#include "Estimation.h"

int	main()
{
	Matrix	A;
	Matrix	R;
	Matrix	A_cpy;	//copy of A matrix for inverse Gauss
	Matrix	T(R.size, 1);
	double	*b_ = NULL;
	double	*x = NULL;

	std::cout << "matrix A:" << std::endl;
	A.print();
	std::cout << "vector b:" << std::endl;
	for (int j = 0; j < A.size; j++)
		std::cout << std::setw(8) << A.rvalue[j] << std::endl;
	std::cout << std::endl;

	//QR
	T.QR_find_x(R);
	if (degenerate_matrix(R) == false)
	{
		b_ = multi_vect(R.rvalue, T);
		x = find_x(R, b_);
		std::cout << "Vector x in QR:" << std::endl;
		for (int j = 0; j < A.size; j++)
		{
			std::cout << std::setw(8) << x[j] << std::endl;
		}
		std::cout << std::endl;

		//find inverse matrix with QR
		Matrix A_(A.size, 1);
		A_.inverse_matrix(R, T);
		std::cout << "Inverse matrix with QR method A_:" << std::endl;
		A_.print();

		//find conditional number with different norm
		std::cout << "Estimation of condition number with cube norme number >= " << cube_estimate_number(A.rvalue, R, T) << std::endl;
		std::cout << "Conditional number with cube norme = " << cube_norm(A_) * cube_norm(A) << std::endl;
		std::cout << std::endl;
		std::cout << "Estimation of condition number with octahedral norme number >= " << octah_estimate_number(A.rvalue, R, T) << std::endl;
		std::cout << "Conditional number with octahedral norme = " << octah_norm(A_) * octah_norm(A) << std::endl;
		std::cout << std::endl;
		std::cout << "Estimation of condition number with octahedral norme number >= " << sfer_estimate_number(A.rvalue, R, T) << std::endl;
		std::cout << "Conditional number with sferical norme = " << sfer_norm(A_) * sfer_norm(A) << std::endl;
		std::cout << std::endl;

		//find residual vector
		std::cout << "Residual vector with cube norme = " << cube_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.size), A.size) << std::endl;
		std::cout << std::endl;
		std::cout << "Residual vector with octahedral norme = " << octah_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.size), A.size) << std::endl;
		std::cout << std::endl;
		std::cout << "Residual vector with sferical norme = " << sfer_vect_norm(diff_vector(A.rvalue, multi_vect(x, A), A.size), A.size) << std::endl;
		std::cout << std::endl;

	}
	else
		std::cout << "Wrong Matrix" << std::endl;

	//delete variable of QR
	if (x)
		delete[] x;
	if (b_)
		delete[] b_;

	//Gauss 
	GaussRight(A);
	if (degenerate_matrix(A) == false)
	{
		GaussLeft(A);
		//find inverse matrix with Gauss
		A_cpy.print();
		Matrix E(A_cpy.size, 1);
		reversematrixgauss(A_cpy, E);
		std::cout << "Inverse matrix with Gauss method A_:" << std::endl;
		E.print();
	}
	else
		std::cout << "Wrong Matrix" << std::endl;

	system("pause");
	return (0);
}