/*******************************************************
Program 4.  Finds the least-squares trig poly S_n of 
specified degree n for a given set of data (x_j,y_j), 
j=0...2m-1.  It is assumed 2 <= n < m and that x_j are 
the std nodes in [-pi,pi], i.e. x_j = -pi + j*pi/m.

Inputs:  data (x_j,y_j),  j=0...2m-1
         degree n (entered at runtime)

Outputs: coeffs a_0...a_n, b_1...b_{n-1} of S_n
         error E_n associated with S_n

Note 1: This program is incomplete; you'll need to 
code the a,b-coefficients and fitting error E as 
indicated below.

Note 2: For any given problem, you'll need to set the 
value of m and specify the input file with the xy-data,
and then compile and run.
*******************************************************/
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;


/*** Specify name of input data file ***/
const int maxChar = 20;
const char myfile[maxChar]="program4.dat" ;


int main() {

  /*** Input/read (x_j,y_j) data ***/
  int m=50;
  double pi=4.0*atan(1.0); // the number pi
  vector x(2*m), y(2*m);

  ifstream fileRead(myfile);
  for(int j=0; j<2*m; j++) {
    fileRead >> x(j) >> y(j); // read xy-pair 
  }

  /*** Echo input file details ***/
  cout << endl;
  cout << "xy-data read from file: " << myfile << endl;
  cout << endl;
  cout << "Number of data points 2m: " << 2*m << endl;
  cout << endl;

  /*** Compute a,b-coefficients ***/
  int n; 
  cout << "Enter trig poly degree n: " << flush;
  cin >> n;
  if( n<2 || n>= m){
    cerr << "Trig LS: bad degree -- exiting" << endl;
    exit(EXIT_FAILURE);
  }

  vector a(n+1), b(n);
  a = 0 ; b = 0 ; // initialize for sum
  

  double asum = 0, bsum = 0;
  //first compute a(0) and a(n) outside loop
  for (int j=0; j<2*m; j++)
  {
    asum += y(j);
  }
  a(0) = (1/m)*asum;
  asum = 0;
  for (int j=0; j<2*m; j++)
  {
    asum += y(j)*cos(n*x(j));
  }
  a(n) = (1/m)*asum;
  asum = 0;

  //now compute a(k) and b(k) simultaneously
  for (int k=1; k<n; k++)
  {
    for (int j=0; j<2*m; j++)
    {
      asum += y(j)*cos(k*x(j));
      bsum += y(j)*sin(k*x(j));
    }
    a(k) = (1.0/m)*asum;
    b(k) = (1.0/m)*bsum;
    asum = 0;
    bsum = 0;
  }


  /*** Print least-squares coeffs to screen ***/
  cout << endl;
  cout << "Trig LS results for n = " << n << ", m = " << m << endl;
  cout << endl;
  cout << "Cosine coeffs: a_0 ... a_n = " << endl;
  cout << a << endl;
  cout << "Sine coeffs: b_1 ... b_{n-1} = " << endl;
  cout << b << endl;

  /*** Compute least-squares fitting error ***/
  double E = 0;
  double S = 0;

  for (int j=0; j<2*m; j++)
  {
    //first compute S(x_j)
    S += a(0)/2.0;
    S += a(n)*cos(n*x(j))/2.0;
    for (int k=1; k<n; k++)
      S += a(k)*cos(k*x(j)) + b(k)*sin(k*x(j));

    //now compute error term
    E += (y(j) - S)*(y(j) - S);
    S = 0;
  }


  cout << endl;
  cout << "Fitting error: E = " << endl;
  cout << E << endl;
    
  return 0; //terminate main program
}