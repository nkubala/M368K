/*******************************************************
Function to implement the steepest descent method to 
minimize a function g(x).

Inputs:  
  geval    Function to evaluate g(x), scalar
  dgeval   Function to evaluate dg/dx(x), n vector
  xk       Initial guess for x, n vector  
  maxIter  Maximum number of descent iterations
  tol      Tolerance parameter for ||dg/dx||

  ainit    Initial guess for alpha (step size)
  maxSrch  Maximum number of alpha search steps


Outputs: 
  xk       Approx local minimum point of g(x)
  k        Number of descent steps taken


Note 1: This function is incomplete; you'll need to code 
the quadratic interpolation step as indicated below.

Note 2: The functions geval and dgeval are assumed 
to be defined externally (e.g. by calling program).
*******************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std ;

/*** Declare user-defined functions to be used ***/
void geval(vector&, double&) ;  //defined in program8.cpp
void dgeval(vector&, vector&) ; //defined in program8.cpp


int descent(vector& xk, int n, int maxIter, 
                    double tol, int maxSrch, double ainit){
  int k=0, kSrch=0 ;  
  double gk, dgkNorm ;
  double a3, ga3, a0, ga0 ;
  double a1, ga1, a2, ga2;
  double h1, h2, h3;
  vector dgk(n), xa3(n), xa0(n) ;
  vector xa2(n);

  geval(xk,gk) ; //evaluate g(xk)
  dgeval(xk,dgk) ;  //evaluate dg(xk)
  dgkNorm = vecL2Norm(dgk) ; 
  
  while(dgkNorm>=tol && k<maxIter) {
    a3 = ainit ;
    xa3 = xk - scaleVec(a3/dgkNorm,dgk) ;
    geval(xa3,ga3) ;

    kSrch = 0 ;
    while( ga3>=gk && kSrch<maxSrch ){ 
      a3 = 0.5*a3 ;
      xa3 = xk - scaleVec(a3/dgkNorm,dgk) ;
      geval(xa3,ga3) ;
      kSrch++ ;
    }
    if(kSrch>=maxSrch){
      cout << "Descent: max iter exceeded in alpha search" << endl;
    }

    a0 = a3 ; ga0 = ga3 ; xa0 = xa3 ;

    /*** compute quadratic interpolation 
         quantities a0, xa0 and ga0 here ***/

    a2 = a3 / 2.0;
    xa2 = xk - scaleVec(a2/dgkNorm,dgk);
    geval(xa2, ga2);

    h1 = (ga2 - gk) / a2;
    h2 = (ga3 - ga2) / (a3 - a2);
    h3 = (h2 - h1) / a3;

    a0 = 0.5*(a2 - (h1/h3));
    xa0 = xk - scaleVec(a0/dgkNorm, dgk);
    geval(xa0, ga0);



    if(a0>=0 && a0<=a3 && ga0<ga3){
      xk = xa0 ;
      gk = ga0 ;
    }
    else {
      xk = xa3 ;
      gk = ga3 ;
    }

    dgeval(xk,dgk) ;
    dgkNorm = vecL2Norm(dgk) ; 
    k++ ;
    cout << "Descent: |Gradg_k|_2 = " << dgkNorm << endl ;
  }


  if(dgkNorm < tol) {
    cout << "Descent: solution converged" << endl ;
  } 
  else {
    cout << "Descent: max iterations exceeded" << endl;
  }

  return k;
}