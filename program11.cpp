/************************************************************
Program 11.  Uses the central-difference method to find an 
approximate solution of an elliptic BVP in a rectangular 
domain of the form

 P uxx + Q uyy + p ux + q uy + r u = f,  a<=x<=b, c<=y<=d  
 u(a,y) = ga(y),     u(b,y) = gb(y),          c<=y<=d  
 u(x,c) = gc(x),     u(x,d) = gd(x),          a<=x<=b  


Inputs:  
  PDEeval       Function to evaluate P,Q,p,q,r,f
  BCeval        Function to evaluate ga,gb,gc,gd
  a,b,c,d       Domain parameters
  N,M           Number of interior x,y pts (N+2,M+2 total pts) 
  x             Grid point vector: x(i)=a+i*dx, i=0...N+1
  y             Grid point vector: y(j)=c+j*dy, j=0...M+1

Outputs: 
  x             Grid point vector: x(i)=a+i*dx, i=0...N+1
  y             Grid point vector: y(j)=c+j*dy, j=0...M+1
  u             Approx soln: u(i,j)=soln at x(i),y(j), 
                i=0...N+1, j=0...M+1


Note 1: The function file linearcd2D.cpp is incomplete; 
you'll need to code the A-matrix and G-vector as 
indicated in that file.

Note 2: For any given problem, the functions PDEeval 
and BCeval must be changed.

Note 3: For any given problem, the grid parameters a,b,
c,d,N,M must be specified.  

Note 4: Gauss elimination is used to solve the system,
so only moderate values of N,M can be used.

Note 5: To compile this program use the command (all 
on one line)

  c++ -o program11  matrix.cpp  gauss_elim.cpp 
                          linearcd2D.cpp  program11.cpp

Note 6: The program output is written to a file.
************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;


/*** Define output file ***/
const char myfile[20]="program11.out" ;
ofstream prt(myfile) ;


/*** Declare external function ***/
int linearcd2D(int, int, double, double, double, double,
                                      vector&, vector&, matrix&) ;


/*** Define P(x,y), Q(x,y), p(x,y), q(x,y), r(x,y), f(x,y) ***/
void PDEeval(const double& x, const double& y, 
             double& P, double& Q, double& p, double& q,
                                           double& r, double& f){

  // P = 1 ;
  // Q = 1 ;
  // p = 0 ;
  // q = x ;
  // r = 0 ;
  // f = x*y ;

  P = 0.3 + 0.5*y;
  Q = 0.3 + 0.5*y;
  p = -(7 - 5*x);
  q = -5*y ;
  r = 0 ;

  f = -(10*exp(-30*pow(x-0.2,2) - 30*pow(y-0.2,2)));
  // f = -(10*exp(-30*pow(x-0.2,2) - 30*pow(y-0.2,2)) + 8*exp(-30*pow(x-0.2,2)-30*pow(y-0.6,2)));
  // f = -(10*exp(-30*pow(x-0.2,2) - 30*pow(y-0.2,2)) + 8*exp(-30*pow(x-0.8,2)-30*pow(y-0.2,2)));
}

/*** Define ga(y), gb(y), gc(x), gd(x) ***/
void BCeval(const double& x, const double& y,
            double& ga, double& gb, double& gc, double& gd){

  // ga = y ;
  // gb = sqrt(y) ;
  // gc = 0 ;
  // gd = 1 ;
  ga = 0.1;
  gb = 0.1;
  gc = 0.1;
  gd = 0.1;

}

int main() {
  /*** Define problem parameters ***/
  int N=29, M=29, success_flag=0 ;  
  matrix u(N+2,M+2) ;
  vector x(N+2), y(M+2) ; 
  double a=0, b=2, c=0, d=1 ; 
  double dx=(b-a)/(N+1), dy=(d-c)/(M+1) ;


  /*** Construct grid ***/
  for(int i=0; i<=N+1; i++){
    x(i) = a + i*dx ;
  }
  for(int j=0; j<=M+1; j++){
    y(j) = c + j*dy ;
  } 


  /*** Call central-difference method ***/
  success_flag=linearcd2D(N,M,a,b,c,d,x,y,u) ;


  //find maximum of u(x,y)
  double max = -1;
  double xmax, ymax;
  int i,j, xidx, yidx;
  for (i=0; i<N+2; i++)
  {
    for (j=0; j<M+2; j++)
    {
      if (u(i,j) > max)
      {
        max = u(i,j);
        xidx = i;
        yidx = j;
      }
    }
  }

  xmax = x(xidx);
  ymax = y(yidx);

  cout << endl;
  cout << "Max is " << max << " at x = " << xmax << ", y = " << ymax << endl;
  cout << endl;


  /*** Print results to output file ***/
  prt.setf(ios::fixed) ;
  prt << setprecision(5) ;
  cout << "Linear-CD-2D: output written to " << myfile << endl ;
  prt << "Linear-CD-2D results" << endl ;
  prt << "Number of interior x-grid pts: N = " << N << endl ;
  prt << "Number of interior y-grid pts: M = " << M << endl ;
  prt << "Approximate solution: x_i, y_j, u_ij" << endl ;
  for(int i=0; i<=N+1; i++){
    for(int j=0; j<=M+1; j++){
      prt << setw(8) << x(i) ;
      prt << "   " ;
      prt << setw(8) << y(j) ;
      prt << "   " ;
      prt << setw(8) << u(i,j) ;
      prt << endl;
    }
  }

  return 0 ; //terminate main program
}
