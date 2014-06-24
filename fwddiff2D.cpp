/************************************************************
Function to implement the forward-difference method to find 
an approximate solution of a parabolic IBVP in a rectangular 
domain of the form 

 ut = P uxx + Q uyy + p ux + q uy + r u + eta, 
                                   a<=x<=b, c<=y<=d, 0<=t<=T  

 u(a,y,t) = ga(y,t),   u(b,y,t) = gb(y,t),  c<=y<=d, 0<=t<=T  
 u(x,c,t) = gc(x,t),   u(x,d,t) = gd(x,t),  a<=x<=b, 0<=t<=T  

 u(x,y,0) = f(x,y),                         a<=x<=b, c<=y<=d  


Inputs:  
  PDEeval       Function to evaluate P,Q,p,q,r,eta
  BCeval        Function to evaluate ga,gb,gc,gd
  a,b,c,d       Space domain parameters
  N,M           Number of interior x,y pts (N+2,M+2 total pts) 
  x             Grid point vector: x(i)=a+i*dx, i=0...N+1
  y             Grid point vector: y(j)=c+j*dy, j=0...M+1
  t,dt          Current time t and time step dt
  u             Approx soln at time t: u(i,j) is soln 
                at x(i),y(j), i=0...N+1, j=0...M+1

Outputs: 
  x             Grid point vector: x(i)=a+i*dx, i=0...N+1
  y             Grid point vector: y(j)=c+j*dy, j=0...M+1
  u             Approx soln at time t+dt: u(i,j) is soln 
                at x(i),y(j), i=0...N+1, j=0...M+1


Note 1: This function is incomplete; you'll need to finish
coding the method as indicated below.  

Note 2: The functions PDEeval and BCeval are assumed 
to be defined externally (e.g. by calling program).

Note 3: Rather than use the single-label index l=1...NM
for the interior xy-grid, it is convenient to use the 
double-label indices i=0...N+1, j=0...M+1 for the entire 
xy-grid.
************************************************************/
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;


/*** Ext fun: PDEeval(x,y,t,P,Q,p,q,r,eta) ***/
void PDEeval(const double&, const double&, double&,
                           double&, double&, double&, 
                             double&, double&, double&) ;

/*** Ext fun: BCeval(x,y,t,gLeft,gRight,gBottom,gTop) ***/
void BCeval(const double&, const double&, double&,
                     double&, double&, double&, double&) ;


/*** Main function: forward-difference method ***/
int fwddiff2D(int N, int M, 
              double a, double b, double c, double d, 
              vector& x, vector& y, double& t, double& dt, 
                                                  matrix& u){

  int success_flag=0 ;
  double P, Q, p, q, r, eta, tn, tnn ;
  double gLeft, gRight, gBottom, gTop ;
  double uLeft=0, uRight=0, uBottom=0, uTop=0, uCenter=0 ;
  double dx=(b-a)/(N+1), dy=(d-c)/(M+1) ; //grid params
  double ahatij=0,bhatij=0,chatij=0,dhatij=0,ehatij=0 ; //fwd-diff params
  matrix unew(N+2,M+2) ; //temporary variable


  /*** Compute u at t_{n+1} at interior grid points ***/
  for(int i=1; i<=N; i++){ //actual i-index on grid
    for(int j=1; j<=M; j++){ //actual j-index on grid

      tn = t ;
      PDEeval(x(i),y(j),tn,P,Q,p,q,r,eta) ;
      BCeval(x(i),y(j),tn,gLeft,gRight,gBottom,gTop) ;
 
      /*** Build the ahat_{i,j},...,dhat_{i,j} params here 
      and compute u^{n+1} in terms of u^{n}, eta-values
      and g-values.  A few lines are given as an example.  
      We use the more description notation uCenter in 
      place of u_{i,j}, uTop in place of u_{i,j+1}, etc. ***/

      ahatij = P/(dx*dx) - p/(2.0*dx) ;
      bhatij = r - 2.0*P/(dx*dx) - 2.0*Q/(dy*dy) ;
      chatij = P/(dx*dx) + p/(2.0*dx) ;
      dhatij = Q/(dy*dy) + q/(2.0*dy);
      ehatij = Q/(dy*dy) - q/(2.0*dy);

      uCenter = u(i,j);

      if( i>1 )
        uLeft = u(i-1,j); //interior value
      else
        uLeft = gLeft; //boundary value

      if( i<N )
        uRight = u(i+1,j); //interior value
      else
        uRight = gRight; //boundary value

      if( j>1 )
        uBottom = u(i,j-1);
      else
        uBottom = gBottom;
      
      if( j<M )
        uTop = u(i,j+1); 
      else
        uTop = gTop;


      unew(i,j) = u(i,j) + dt*dhatij*uTop
                         + dt*ahatij*uLeft 
                         + dt*bhatij*uCenter 
                         + dt*chatij*uRight
                         + dt*ehatij*uBottom
                         + dt*eta ; //fwd-diff formula
    }
  }
 
  /*** Compute u at t_{n+1} at boundary grid points ***/
  for(int i=0; i<=N+1; i++){ //actual i-index on grid
    tnn = t + dt ;
    BCeval(x(i),y(M+1),tnn,gLeft,gRight,gBottom,gTop) ;
    unew(i,M+1) = gTop ; 
    BCeval(x(i),y(0),tnn,gLeft,gRight,gBottom,gTop) ;
    unew(i,0) = gBottom ; 
  }
  for(int j=0; j<=M+1; j++){ //actual j-index on grid
    tnn = t + dt ;
    BCeval(x(N+1),y(j),tnn,gLeft,gRight,gBottom,gTop) ;
    unew(N+1,j) = gRight ; 
    BCeval(x(0),y(j),tnn,gLeft,gRight,gBottom,gTop) ;
    unew(0,j) = gLeft ; 
  }

  u = unew ; //update u-matrix with new values

  return success_flag=1 ;
}
