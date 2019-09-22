#include <iostream>
#include <fstream>
#include <cmath>

class Matrix
{
public:
  double *value;
  double *rvalue;
  double **mymap;
  int size;
public:
  Matrix ()
  {
    int i = 1;
    value = new double[size*size];
    mymap = new double*[size];
    rvalue = new double[size];
    while (*mymap)
    {
      *mymap = &value[size*i];
      mymap++;
    }
    double number;
    std::ifstream file;
    file.open("matrixexample");
    while (!file.eof())
    {
      for (int i = 0; i != size; i++)
      {
        file >> number;
        *value = number;
        value++;
      }
      file >> number;
      *rvalue = number;
      rvalue++;
    }
    file.close();
  }

  Matrix (Matrix& A)
  {
    Matrix C;
    int i = 1;
    C.value = new double[A.size*A.size];
    C.mymap = new double*[A.size];
    C.rvalue = new double[size];
    while (*value)
    {
      *C.value=*(A.value);
      C.value++;
      (A.value)++;
    }
    while (*C.mymap)
    {
      *C.mymap = &C.value[size*i];
      C.mymap++;
      C.rvalue = A.rvalue;
      A.rvalue++;
      C.rvalue++;
    }
    return C;
  }

  void Tranc(Matrix& A)
  { // i - string, j - column
    double temp;
    for (int i = 0; i != A.size; i++)
    {
      for (int j =0; j != A.size; j++)
      {
        if (j>i)
        {
          temp = A.value[(A.mymap[i]-A.value)+j];
          A.value[(A.mymap[i]-A.value)+j] = A.value[(A.mymap[i]-A.value)+j + i - 1];
          A.value[(A.mymap[i]-A.value)+j + i - 1] = temp;
        }
      }
    }
  }

Matrix operator * (Matrix A, Matrix B)
{

  ~Matrix()
  {
    delete[] value;
    delete[] mymap;
    delete[] rvalue;
  }
};
