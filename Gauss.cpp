#include <iostream>
#include <cmath>
#include "MatrixClass.h"

// search max and change strings
Matrix delim(Matrix &A)
{
  for (int i = 0; i != (A.size * A.size); i++)
  {
    A.value[i]=A.value[i]/A.value[A.mymap[i / A.size] - A.value];
    A.rvalue[i]=A.rvalue[i]/A.value[A.mymap[i / A.size] - A.value];
  }
  return A;
}

Matrix vych(Matrix A, int k)
{
  //strings
  for (int i = k+1; i != A.size; i++)
  {
    //columns
    for (int j = 0; j != A.size; j++)
    {
      A.value[(A.mymap[i] - A.value)+j]-= A.value[(A.mymap[k] - A.value)+j];
    }
    A.rvalue[i] -= A.rvalue[k];
  }
  return A;
}

void GaussRight(Matrix& A)
{
  // for k-string our matrix
  for (int k = 0; k != A.size; k++)
  {
    std::cout << "!\n";
    //search max in column
    double max = A.value[0];
    int maxstring = 0;
    for (int i = 0; i != A.size; i++)
    {
      if (A.value[i] > max)
      {
        max = A.value[i];
        maxstring = i / A.size;
      }
    }
    std::cout << "!!!!!\n";
    A.mymap[k] = &A.value[maxstring*A.size];
    A.mymap[maxstring*A.size]= &A.value[k*A.size];
    A = delim(A);
    A = vych(A, k);
  }
  std::cout <<"GaussRight right!\n";
}

void GaussLeft(Matrix A, double *x)
{
  for (int i = A.size -1; i >= 0; i--)
  {
    double sum = 0;
   for (int j= i+1;j < A.size; j++)
   {
     sum += A.value[(A.mymap[i]-A.value)+j]*x[j];
   } 
   x[i] = (A.rvalue[i] - sum)/A.value[(A.mymap[i]-A.value)+i-1]; 
  }
}

int main()
{
  Matrix A;
  GaussRight(A);
  double *x = new double[4];
  GaussLeft(A, x);
  for (int i = 0; i < 4; i++)
  std::cout<< x[i] << "  ";
  delete[] x;

}
