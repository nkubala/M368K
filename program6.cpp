/*******************************************************
Program 6.  Uses the QR method to find the eigenvalues 
and eigenvectors of a symmetric, tri-diagonal matrix A.

Inputs:  
   A        nxn symmetric, tri-diagonal matrix
   maxIter  iteration control param
   tol      tolerance parameter

Outputs: 
   D        nxn diagonal matrix of eigenvalues
   Q        nxn orthogonal matrix of eigenvectors

Note 1: The function file qr.cpp is incomplete; 
you'll need to code the P-matrix and its transpose 
as indicated in that file.

Note 2:  For any given problem, you'll need to 
specify the values of {n,maxIter,tol} and input
the matrix A.
*******************************************************/
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;


/*** Declare user-defined functions to be used ***/
int qr(matrix&, matrix&, matrix&, int, int, double) ;


int main() {
  /*** Define and input problem data ***/
  int n=4, maxIter=50, iter ;  
  double tol=1e-5 ;
  matrix A(n,n), D(n,n), Q(n,n) ;

  // A(0,0)= 3; A(0,1)=-1; A(0,2)= 0; A(0,3)= 0; 
  // A(1,0)=-1; A(1,1)= 5; A(1,2)=0.5; A(1,3)= 0;
  // A(2,0)= 0; A(2,1)=0.5; A(2,2)= 4; A(2,3)=-1; 
  // A(3,0)= 0; A(3,1)= 0; A(3,2)=-1; A(3,3)= 2; 
  A(0,0)=30; A(0,1)=-14.142; A(0,2)=0; A(0,3)=0;
  A(1,0)=-14.142; A(1,1)=17.5; A(1,2)=-6.124; A(1,3)=0;
  A(2,0)=0; A(2,1)=-6.124; A(2,2)=10.002; A(2,3)=-3.873;
  A(3,0)=0; A(3,1)=0; A(3,2)=-3.873; A(3,3)=3;

  /*** Print data to screen ***/
  cout << endl ; 
  cout << "Given: A = " << endl ;
  cout << A << endl ;

  /*** Set a fixed precision for printing output ***/
  cout.setf(ios::fixed) ;
  cout << setprecision(6) ;

  /*** Call the QR function ***/
  iter = qr(A,Q,D,n,maxIter,tol) ;

  /*** Print results to screen ***/
  cout << endl;
  cout << "Iteration index: k = " << iter << endl;
  cout << "Approximate eigval matrix: D^(k) = " << endl;
  cout << D << endl ;
  cout << "Approximate eigvec matrix: Q^(k) = " << endl;
  cout << Q << endl ;

  return 0 ; //terminate main program
}