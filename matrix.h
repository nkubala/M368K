/*****************************************************************
General purpose library: header file

General purpose library of vector and matrix data types and 
related functions.   See file matrix.cpp for implementation.

Usage notes:

  vector v(n)   defines a vector v of n doubles, 
                indexed from 0...n-1.

  matrix A(m,n) defines a matrix A of mxn doubles,
                indexed from 0...m-1 and 0...n-1.

  dim(v)        returns dimension of vector v
  dim(A,0)      returns dimension m of matrix A
  dim(A,1)      returns dimension n of matrix A
  v=w           sets vector v equal to vector w
  A=B           sets matrix A equal to matrix B
  v=c           sets all entries of vector v equal to scalar c
  A=c           sets all entries of matrix A equal to scalar c
  v+w           adds vectors v and w
  A+B           adds matrices A and B
  v-w           subtracts vectors v and w
  A-B           subtracts matrices A and B
  scaleVec(c,v) returns vector cv, where c scalar, v vector
  scaleMat(c,A) returns matrix cA, where c scalar, A matrix
  matMatMult(A,B) returns matrix AB, where A matrix, B matrix
  matVecMult(A,v) returns vector Av, where A matrix, v vector 
  vecL2Norm(v)  returns L_2 norm of vector v
  vecMaxNorm(v) returns L_infty (max) norm of vector v
  matMaxNorm(A) returns L_infty (max) norm of matrix A
  vecDot(v,w)   returns dot product of vectors v, w
  cout << v     prints entries of vector v to screen
  cout << A     prints entries of matrix A to screen
  cin >> v      loads entries into vector v from keyboard
  cin >> A      loads entries into matrix A from keyboard
  
******************************************************************/

#ifndef MATRIX_H  //compile library only if 
#define MATRIX_H  //it is not defined elsewhere

#include <iostream>
using namespace std;

/**************************************** 
matrix class
*****************************************/
class matrix {
 private:
  int size[2]; //size of matrix
  double* array; //pointer to elements of matrix

 public:
  matrix() ; //default constructor 
  matrix(int,int) ; //constructor 
  ~matrix(); //destructor 
  matrix(const matrix&) ; //copy constructor 
  double operator() (int,int) const ; //get-value index operator 
  double& operator() (int,int) ; //set-value index operator 
  matrix& operator=(const matrix&) ; //equate-two-mats operator
  matrix& operator=(double) ;        //fill-mat-by-scalar operator
  matrix operator+(const matrix&) ;  //add-two-mats operator
  matrix operator-(const matrix&) ;  //subtract-two-mats operator
  friend int dim(const matrix&,int) ; //mat dimension function
  friend matrix scaleMat(const double, const matrix&); //scale mat function 
  friend matrix matMatMult(const matrix&, const matrix&); //mat*mat function
};

/**************************************** 
vector class
*****************************************/
class vector {
 private:
  int size; //size of vector
  double* array; //pointer to elements of vector

 public:
  vector() ; //default constructor 
  vector(int) ; //constructor 
  ~vector(); //destructor 
  vector(const vector&) ; //copy constructor 
  double operator() (int) const ; //get-value index operator 
  double& operator() (int) ; //set-value index operator 
  vector& operator=(const vector&) ; //equate-two-vecs operator
  vector& operator=(double) ;        //fill-vec-by-scalar operator
  vector operator+(const vector&) ;  //add-two-vecs operator
  vector operator-(const vector&) ;  //subtract-two-vecs operator
  friend int dim(const vector&) ;    //vec dimension function
  friend vector scaleVec(const double, const vector&); //scale vec function 
  friend vector matVecMult(const matrix&, const vector&); //mat*vec function
};

/*********************************************** 
input, output operators for vectors and matrices
***********************************************/
ostream& operator<< (ostream&, const vector&); // output a vector
ostream& operator<< (ostream&, const matrix&); // output a matrix by rows
istream& operator>> (istream&, vector&);       // input a vector
istream& operator>> (istream&, matrix&);       // input a matrix by rows

/*********************************************** 
utility functions for vectors and matrices
***********************************************/
double vecL2Norm(const vector&);  // l2-norm of a vector
double vecMaxNorm(const vector&); // l-infinity (max) norm of a vector
double matMaxNorm(const matrix&); // l-infinity (max) norm of a matrix
double vecDot(const vector&, const vector&); // vector dot product

#endif