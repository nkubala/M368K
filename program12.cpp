/************************************************************
Program 12.  Uses the forward-difference method to find an 
approximate solution of a parabolic IBVP in a rectangular 
domain of the form

 ut = P uxx + Q uyy + p ux + q uy + r u + eta, 
                                   a<=x<=b, c<=y<=d, 0<=t<=T  

 u(a,y,t) = ga(y,t),   u(b,y,t) = gb(y,t),  c<=y<=d, 0<=t<=T  
 u(x,c,t) = gc(x,t),   u(x,d,t) = gd(x,t),  a<=x<=b, 0<=t<=T  

 u(x,y,0) = f(x,y),                         a<=x<=b, c<=y<=d  


Inputs:  
  PDEeval       Function to evaluate P,Q,p,q,r,eta
  BCeval        Function to evaluate ga,gb,gc,gd
  ICeval        Function to evaluate f
  a,b,c,d       Space domain parameters
  N,M           Number of interior x,y pts (N+2,M+2 total pts) 
  x             Grid point vector: x(i)=a+i*dx, i=0...N+1
  y             Grid point vector: y(j)=c+j*dy, j=0...M+1
  T             Time interval parameter (final time) 
  L             Number of time steps

Outputs: 
  x             Grid point vector: x(i)=a+i*dx, i=0...N+1
  y             Grid point vector: y(j)=c+j*dy, j=0...M+1
  u             Approx soln at time t=T: u(i,j) is soln 
                at x(i),y(j), i=0...N+1, j=0...M+1


Note 1: The function file fwddiff2D.cpp is incomplete; you'll 
need to finish coding the method as indicated in that file.

Note 2: For any given problem, the functions PDEeval, BCeval 
and ICeval must be changed.

Note 3: For any given problem, the grid parameters a,b,c,d,T
and N,M,L must be specified.  

Note 4: To compile this program use the command (all on one 
line)

  c++ -o program12  matrix.cpp  fwddiff2D.cpp  program12.cpp

Note 5: The program output is written to a file.
************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;


/*** Define output file ***/
const char myfile[20]="program12.out" ;
ofstream prt(myfile) ;


/*** Declare external function ***/
int fwddiff2D(int, int, double, double, double, double,
                   vector&, vector&, double&, double&, matrix&) ;

/*** Define P(x,y,t), Q(x,y,t), 
                     p(x,y,t), q(x,y,t), r(x,y,t), eta(x,y,t) ***/
void PDEeval(const double& x, const double& y, double& t,
                    double& P, double& Q, double& p, double& q,
                                          double& r, double& eta){


  P = 0.03; 
  Q = 0.03; 
  p = -y; 
  q = x;
  r = 0 ; 
  eta = 0.6 * exp(-30*pow(x-0.2,2) - 30*pow(y-0.2,2))
      - 0.7 * exp(-30*pow(x-0.6,2) - 30*pow(y-0.2,2))
      - 0.6 * exp(-30*pow(x+0.5,2) - 30*pow(y,2));

}

/*** Define ga(y,t), gb(y,t), gc(x,t), gd(x,t) ***/
void BCeval(const double& x, const double& y, double& t,
            double& ga, double& gb, double& gc, double& gd){


  ga = 1 ; 
  gb = 1 ;
  gc = 1 ; 
  gd = 1 ;  

}

/*** Define f(x,y) ***/
void ICeval(const double& x, const double& y, double& f){

  f = 1 ;
}


int main() {
  /*** Define problem parameters ***/
  //int N=19, M=19, L=100, success_flag=0 ;  //set 1
  //int N=29, M=29, L=125, success_flag=0 ;  //set 2
  //int N=39, M=39, L=250, success_flag=0 ;  //set 3

  int N=29, M=29, success_flag=0 ;  //part c

  matrix u(N+2,M+2) ;
  vector x(N+2), y(M+2) ; 
  double a=-1, b=1, c=-1, d=1 ; 
  double dx=(b-a)/(N+1), dy=(d-c)/(M+1) ;

  // double T=20, dt=T/L ; 
  // int L=50; double T=0.5, dt=T/L;
  // int L=300; double T=3, dt=T/L;
  // int L=1000; double T=10, dt=T/L;
  // int L=2000; double T=20, dt=T/L;
  int L=4000; double T=40, dt=T/L;


  /*** Construct xy-grid ***/
  for(int i=0; i<=N+1; i++){
    x(i) = a + i*dx ;
  }
  for(int j=0; j<=M+1; j++){
    y(j) = c + j*dy ;
  } 

  /*** Load initial condition ***/
  double t=0 ;
  double gLeft, gRight, gBottom, gTop, f ;
  for(int i=0; i<=N+1; i++){ //actual i-index on grid
    for(int j=0; j<=M+1; j++){ //actual j-index on grid
      BCeval(x(i),y(j),t,gLeft,gRight,gBottom,gTop) ;
      ICeval(x(i),y(j),f) ;
      if( j==M+1 ){ u(i,j) = gTop ; }
      if( i==N+1 ){ u(i,j) = gRight ; }
      if( i==0 ){ u(i,j) = gLeft ; }
      if( j==0 ){ u(i,j) = gBottom ; }
      if((i>0)&&(i<N+1)&&(j>0)&&(j<M+1)){ u(i,j) = f ; }
    }
  }

  /*** Call fwd-diff method at each time step (overwrites u) ***/
  for(int n=0; n<L; n++){
    success_flag = fwddiff2D(N,M,a,b,c,d,x,y,t,dt,u) ;
    t = t + dt ;
  }

  //compute uavg
  double uavg = 0;
  for (int i=0; i<N+1; i++)
  {
    for (int j=0; j<M+1; j++)
    {
      uavg += u(i,j);
    }
  }
  uavg /= (double)((N+1)*(M+1));
  printf("u average is %f\n", uavg);

  /*** Print results at final time to output file ***/
  prt.setf(ios::fixed) ;
  prt << setprecision(5) ;
  cout << "Fwd-Diff-2D: output written to " << myfile << endl ;
  prt << "Fwd-Diff-2D results" << endl ;
  prt << "Number of interior x-grid pts: N = " << N << endl ;
  prt << "Number of interior y-grid pts: M = " << M << endl ;
  prt << "Number of time steps: L = " << L << endl ;
  prt << "Final time of simulation: t = " << t << endl ;
  prt << "Average concentration: u_avg = " << uavg << endl;
  prt << "Approximate solution at time t: x_i, y_j, u_ij" << endl ;
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
