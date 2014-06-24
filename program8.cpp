/*******************************************************
Program 8.  Uses the steepest descent method to find a 
local minimum of a function g(x).

Inputs:  
  geval    Function to evaluate g(x), scalar
  dgeval   Function to evaluate dg/dx(x), n vector
  x        Initial guess for x, n vector  
  maxIter  Maximum number of descent iterations
  tol      Tolerance parameter for ||dg/dx||

  ainit    Initial guess for alpha (step size)
  maxSrch  Maximum number of alpha search steps


Outputs: 
  x        Approx local minimum point of g(x)
  g        Approx value of g(x) at local min
  k        Number of descent steps taken


Note 1: For any given problem the functions geval and
dgeval must be changed.  These functions compute the
scalar g and vector dg for any given x.

Note 2: The function file descent.cpp is incomplete; 
you'll need to code the quadratic interpolation step
as indicated in that file.

Note 3: To compile this program use the command (all 
on one line)

  c++ -o program8 matrix.cpp descent.cpp program8.cpp

*******************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std ;


/*** Declare function with steepest descent algorithm ***/
int descent(vector&, int, int, double, int, double) ;


/*** Define g function for problem ***/
void geval(vector& x, double& g){
  //double a=1, b=-0.1 ; //define any constants
  double a=1.0, b=1.5, l=0.7;
  double pi=4.0*atan(1.0) ; //the number pi
  double y, z;

  // g = pow(a,2) / (pow(l,4) + pow(x(0),4) - 4*pow(l,3)*x(0)*sin(x(1))
  //  - 4*l*pow(x(0),3)*sin(x(1)) + 4*pow(l,2)*pow(x(0),2)*pow(sin(x(1)),2));
  // g -= (2*a) / (pow(l,2) + pow(x(0),2) - 2*l*x(0)*sin(x(1)));
  // g += pow(b,2) / (pow(l,4) + pow(x(0),4) + 4*pow(l,3)*x(0)*sin(x(1))
  //  + 4*l*pow(x(0),3)*sin(x(1)) + 4*pow(l,2)*pow(x(0),2)*pow(sin(x(1)),2));
  // g -= (2*b) / (pow(l,2) + pow(x(0),2) + 2*l*x(0)*sin(x(1)));

  //cleaner to evaluate y and z first
  y = pow((l*cos(x(1))),2) + pow((x(0)-l*sin(x(1))),2);
  z = pow((l*cos(x(1))),2) + pow((x(0)+l*sin(x(1))),2);

  g = (pow(a,2) / pow(y,2)) - ((2*a) / y) + (pow(b,2) / pow(z,2)) - ((2*b) / z);
}


/*** Define dg function (gradient) for problem ***/
void dgeval(vector& x, vector& dg){
  // double a=1, b=-0.1 ; //define any constants
  double a=1.0, b=1.5, l=0.7;
  double pi=4.0*atan(1.0) ; //the number pi
  double r = x(0), s = sin(x(1)), c = cos(x(1));

  // dg(0) = b*x(1) + 2*x(0) ;
  // dg(1) = b*x(0) + 2*x(1) ;

  //this is horrible.
  dg(0) = 
    (pow(a,2) * ((-4*pow(r,-5) - pow(l,-2)*pow(r,-3) + (1/4)*pow(l,-3)*pow(r,-2)*pow(s,-1)
    + (3/4)*pow(l,-1)*pow(r,-4)*pow(s,-1) - (1/2)*pow(l,-2)*pow(r,-3)*pow(s,-2))))

    - (2*a)*(-2*pow(r,-3) + (1/2)*pow(r,-2)*pow(l,-1)*pow(s,-1))
    - (2*b)*(-2*pow(r,-3) - (1/2)*pow(r,-2)*pow(l,-1)*pow(s,-1))

    + (pow(b,2) * ((-4*pow(r,-5) - pow(l,-2)*pow(r,-3) - (1/4)*pow(l,-3)*pow(r,-2)*pow(s,-1)
    - (3/4)*pow(l,-1)*pow(r,-4)*pow(s,-1) - (1/2)*pow(l,-2)*pow(r,-3)*pow(s,-2))));

  dg(1) = 
    (pow(a,2) * ((1/4)*pow(l,-3)*pow(r,-1)*c*pow(s,-2) + (1/4)*pow(l,-1)*pow(r,-3)*c*pow(s,-2) - (1/2)*pow(l,-2)*pow(r,-2)*c*pow(s,-3)))

    - (2*a)*((1/2)*pow(l,-1)*pow(r,-1)*c*pow(s,-2))
    - (2*b)*((-1/2)*pow(l,-1)*pow(r,-1)*c*pow(s,-2))

    + (pow(b,2) * ((-1/4)*pow(l,-3)*pow(r,-1)*c*pow(s,-2) + (-1/4)*pow(l,-1)*pow(r,-3)*c*pow(s,-2) - (1/2)*pow(l,-2)*pow(r,-2)*c*pow(s,-3)));
}


int main() {
  /*** Define problem parameters ***/
  int n=2, maxIter=100, iter=0, maxSrch=20 ;
  double tol=1e-6, ainit=0.2, g ;
  vector x(n) ;
  x(0) = 1.0 ; x(1) = 0.3 ; //initial guess

  /*** Print problem info ***/
  cout << setprecision(10) ;
  cout << endl ;
  cout << "System size: n = " << n << endl ;
  cout << "Initial guess: x^(0) = " << endl ;
  cout << x << endl ;

  /*** Call steepest descent w/initial guess***/
  iter = descent(x,n,maxIter,tol,maxSrch,ainit) ;

  /*** Print results to screen ***/
  geval(x,g) ;
  cout << endl ;
  cout << "Iteration index: k = " << iter << endl ;
  cout << endl ;
  cout << "Approx solution: x = " << endl ;
  cout << x << endl ;
  cout << "Approx g-value: g = " << g << endl ;

  return 0; //terminate main program
}