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
    std::ifstream file;
    file.open("matrixexample.txt");
    file >> size;
    int i = 0;
    value = new double[size*size];
    mymap = new double*[size];
    rvalue = new double[size];
    while (*mymap)
    {
      *mymap = &value[size*i];
      mymap++;
    }
    double number;
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
    std::cout<< "Class ready to go\n";
    file.close();
  }

 /* Matrix (Matrix& A)
{
    Matrix C;
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
      *C.mymap = &C.value[size];
      C.mymap++;
      C.rvalue = A.rvalue;
      A.rvalue++;
      C.rvalue++;
    }
  } */

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
  void print(Matrix A)
  {
    for (int i = 0; i != (A.size*A.size); i++)
    {
      std::cout << A.value[i];
      if (i % (A.size - 1) == 0)
      {
        std::cout << "\n";
      }
    }
  }



  ~Matrix()
  {
    delete[] value;
    delete[] mymap;
    delete[] rvalue;
  }
};
