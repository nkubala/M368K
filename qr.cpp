/*******************************************************************
Function to find all the eigenvalues and eigenvectors of a 
symmetric, tri-diagonal matrix A using the QR method with 
zero shifts.  Specifically, an orthogonal matrix Q=[v1,...,vn] 
and a diagonal matrix D=diag(lam1,...,lamn) are found such that, 
approximately,
                            Q^T A Q = D.


Inputs: matrix A(n,n), matrix Q(n,n), matrix D(n,n),
        integer n, integer maxIter, double tol.

Outputs: matrix Q(n,n), matrix D(n,n), integer iter (returned)


Note 1: This function is incomplete; you'll need to code the 
P-matrix and its transpose as indicated below.

Note 2: The matrix A is not over-written by this function.

Note 3: The matrices Q and D are over-written; their initial
values are not used; they are re-initialized below.
********************************************************************/
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;


int qr(matrix& A, matrix& Q, matrix& D, 
                 int n, int maxIter, double tol) {
  
  int iter=0 ;
  double error=1, c, s, theta; //for calculating rotation matrix
  double pi=4.0*atan(1.0) ; //the number pi
  matrix B(n,n), P(n,n), Pt(n,n), G(n,n), R(n,n), Eye(n,n) ;

  Eye = 0 ;
  for(int i=0; i<n; i++) {
    Eye(i,i) = 1 ; //the identity matrix of size n
  }

  B = A ; //work with copy of A to preserve original matrix
  D = 0 ;
  Q = Eye ;
  
  while(iter<maxIter && error>=tol) {
    G = Eye ;
    for(int i=1; i<n; i++){
      P  = Eye ; 
      Pt = Eye ;

      //calculate theta, and cos/sin
      theta = atan(B(i,i-1)/B(i-1,i-1));
      c = cos(theta);
      s = sin(theta);

      //set values in rotation matrices
      P(i-1,i-1)=c; P(i-1,i)=s;
      P(i,i-1)=-s; P(i,i)=c;

      Pt(i-1,i-1)=c; Pt(i-1,i)=-s;
      Pt(i,i-1)=s; Pt(i,i)=c;
      

      B = matMatMult(P,B) ; //create zero at {i,i-1} 
      G = matMatMult(G,Pt) ; //accumulate Pt 
    }

    B = matMatMult(B,G) ; //finish similarity transformation
    Q = matMatMult(Q,G) ; //update current eigenvector matrix
    for(int i=0; i<n; i++) {
      D(i,i) = B(i,i) ; //extract current diagonal matrix
    }
    R = B - D ; //form residual of off-diagonal portion
    error = matMaxNorm(R) ;
    iter++ ;
  }

  if(error < tol) {
    cout << "QR method: soln converged, |R|_inf = " << error << endl;
  }
  else {
    cout << "QR method: max iterations exceeded" << endl;
  }
  return iter;
}