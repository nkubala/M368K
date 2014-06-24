#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;

int genpower(matrix&, vector&, double&, int, int, double) ;

int main()
{
	matrix A(3,3);
	A(0,0)=1; A(0,1)=0.7103; A(0,2)=-3.3762;
	A(1,0)=1; A(1,1)=1.7103; A(1,2)=-3.3762;
	A(2,0)=0; A(2,1)=-0.2897; A(2,2)=0.3762;

	vector x(3);
	x=2;


	int maxIter=35, n=3, iter;  
  	double tol=1e-4, lambda; 

	iter=genpower(A,x,lambda,n,maxIter,tol);

	/*** Print results to screen ***/
	cout << "Iteration index: k = " << iter << endl;
	cout << "Approximate eigval: lambda^(k) = " << lambda << endl;
	cout << "Approximate eigvec: x^(k) = " << endl;
	cout << x << endl ;

	return 0; //terminate main program
}