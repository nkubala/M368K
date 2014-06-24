/*******************************************************
Program 5. Finds the mean vector mu in R^n and 
covariance matrix A in R^nxn for a given set of 
data points p^(s) in R^n, s=1...m.  The dominant 
eigenpair (lambda,x) of A are then approximated 
using a power method.

Inputs for finding mu and A: 
p^(s), s=1...m

Inputs for finding largest eigenpair (lambda,x) of A:  
A, x^(0), maxIter, tol

Outputs of power method: 
lambda^(k), x^(k), iteration count k

Note 1: This program is incomplete; you'll need to 
write code to build the mean vector mu and covariance 
matrix A as indicated below.

Note 2:  For any given problem, you'll need to specify
the input data file, and set the values of {n,m,x^(0),
maxIter,tol}.
*******************************************************/
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;


/*** Specify name of input data file ***/
const int maxChar = 20;
const char myfile[maxChar]="program5.dat" ;


/*** Declare user-defined functions to be used ***/
int genpower(matrix&, vector&, double&, int, int, double) ;
int sympower(matrix&, vector&, double&, int, int, double) ;


int main() {
 /*** Read data points from file ***/
  int n=4, m=100;
  matrix datapt(n,m); //  matrix for p_i^(s) points

  ifstream fileRead(myfile);
  for(int s=0; s<m; s++) { // read each row for s=0....m-1
    for (int i=0; i<n; i++) { // read along row for i=0...n-1
      fileRead >> datapt(i,s) ; 
    }
  }

  /*** Echo input file details ***/
  cout << endl;
  cout << "***Input data details***" << endl;
  cout << endl;
  cout << "Data file: " << myfile << endl;
  cout << "Number of data points: m = " << m << endl;
  cout << endl;


  /*** Build mean mu and covariance A ***/
  vector mu(n) ;
  matrix A(n,n) ;
  mu = 0 ; A = 0 ; // initialize for sum


  for (int i=0; i<n; i++) //each entry mu_i
  {
    for (int j=0; j<m; j++) //sum over m data points
    {
      mu(i) += datapt(i,j);
    }
    mu(i) /= m;
  }

  for (int i=0; i<n; i++)
  {
    for (int j=0; j<n; j++)
    {
      for (int k=0; k<m; k++)
      {
        A(i,j) += ( (datapt(i,k) - mu(i)) * (datapt(j,k) - mu(j)) );
      }
      A(i,j) /= m;
    }
  }



  /*** Print mu and A to screen ***/
  cout << "***Results for mu and A***" << endl;
  cout << endl;
  cout << "Mean vector mu = " << endl;
  cout << mu << endl;
  cout << "Covariance matrix A = " << endl;
  cout << A << endl;

  /*** Parameters for power method ***/
  int maxIter=35, iter;  
  double tol=1e-4, lambda; 
  vector x(n);

  x=2; //initial vec
  if(matMaxNorm(A)==0) {A=1;} //reset A if zero

  /*** Print data to screen ***/
  cout << "***Results for power method***" << endl;
  cout << endl; 
  cout << "Given: A = " << endl;
  cout << A << endl;
  cout << "Given: x^(0) = " << endl;
  cout << x << endl;

  /*** Call general power function ***/
  //iter=genpower(A,x,lambda,n,maxIter,tol);
  iter=sympower(A,x,lambda,n,maxIter,tol);
  
  /*** Print results to screen ***/
  cout << "Iteration index: k = " << iter << endl;
  cout << "Approximate eigval: lambda^(k) = " << lambda << endl;
  cout << "Approximate eigvec: x^(k) = " << endl;
  cout << x << endl ;

  return 0; //terminate main program
}