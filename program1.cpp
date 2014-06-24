/****************************************************
Program 1.  Solves Ax=b using Jacobi method.  

Inputs: A, b, x^(0), maxIter, tol
Outputs: x^(k), iteration count k

Here's how to get started:

1) Copy program1.cpp (this file), jacobi.cpp, 
matrix.cpp and matrix.h into your working directory.

2) Compile (and link) the files by typing
"c++ -o program1 matrix.cpp jacobi.cpp program1.cpp"
at the Linux prompt.

/*********************************************/
/*****ADD gauss_seidel.cpp TO FILES TO LINK *****
/*********************************************

3) Type "program1" to run the program.

4) For any given problem, you'll need to set the 
values of the variables {n,A,b,x,maxIter,tol} below, 
and then compile and run as described above.
****************************************************/
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;

/*** Declare user-defined functions to be used ***/
int jacobi(matrix&, vector&, vector&, int, double);

int gauss_seidel(matrix&, vector&, vector&, int, double);


/*** Main program ***/
int main() {

  /*** Define and input problem data ***/
  // int n=3, maxIter=25, iter;  
  // double tol=1e-3;
  // matrix A(n,n);
  // vector x(n), b(n);

  // A(0,0)= 2; A(0,1)=-1; A(0,2)=-1; b(0)=sqrt(5.0);
  // A(1,0)= 1; A(1,1)= 3; A(1,2)=-1; b(1)=exp(0.1);
  // A(2,0)=-1; A(2,1)= 2; A(2,2)= 4; b(2)=2.5;
  // x=0 ; //initial x

  // int n=4, maxIter=50, iter;
  // double tol=1e-3;
  // matrix A(n,n);
  // vector x(n), b(n);

  // A(0,0)=4; A(0,1)=1; A(0,2)=-1; A(0,3)=1; b(0)=-2;
  // A(1,0)=1; A(1,1)=4; A(1,2)=-1; A(1,3)=-1; b(1)=-1;
  // A(2,0)=-1; A(2,1)=-1; A(2,2)=5; A(2,3)=1; b(2)=0;
  // A(3,0)=1; A(3,1)=-1; A(3,2)=1; A(3,3)=3; b(3)=1;
  // x=1000000;


  /***** For Homework Assignment 2 *****/
  //Build input matrices using loop
  int n=100, maxIter=500, iter;
  double tol=1e-5;
  matrix A(n,n);
  vector x(n), b(n);

  int i, j;
  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
    {
      if(j == (i-1))
        A(i,j) = -1;
      else if (j == i)
        A(i,j) = 2.0 + ((double)i / 10.0);
      else if (j == (i+1))
        A(i,j) = -1;
      else
        A(i,j) = 0;
    }
    b(i) = 1 + ((double)i / 20.0);
  }

  /*** Print data to screen ***/
  cout << endl; 
  cout << "Given: A = " << endl;
  cout << A << endl;
  cout << "Given: b = " << endl;
  cout << b << endl;
  cout << "Given: x^(0) = " << endl;
  cout << x << endl;


  /*** Call Jacobi function ***/
  iter=jacobi(A,b,x,maxIter,tol);
  

  /*** Print results to screen ***/
  cout << endl; 
  cout << "Iteration index: k = " << iter << endl;
  cout << endl; 
  cout << "Approximate solution: x^(k) = " << endl;
  cout << x << endl;


  //Call Gauss-Seidel
  x=1000000;
  iter = gauss_seidel(A,b,x,maxIter,tol);

  /*** Print results to screen ***/
  cout << endl; 
  cout << "Iteration index: k = " << iter << endl;
  cout << endl; 
  cout << "Approximate solution: x^(k) = " << endl;
  cout << x << endl;

  return 0;
}
