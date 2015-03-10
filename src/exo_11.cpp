/*__________________________________________________________________________________________________

  autor: Pedro Guarderas
  email: ajusworkopensource@gmail.com
  date: 10-03-2015
  file: exo_11.cpp
  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation;
  either version 2 of the License, or (at your option) any later version.
  __________________________________________________________________________________________________
*/

#include <cmath>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>

#include "linear.hpp"
#include "stosolver.hpp"


using namespace std;


double u( double x, double t ) {
  return sin( 3.14159 * x );
}

double s( double x, double t ) {
  return 0.5*sqrt(x);
}

double psi( double x ) {
  double K = 1.0;
  return max( 0.0, x-K );
}

double V( double x, double t ) {
  return 1.1;
}

double g( Vector< double > X, Vector< double > t ) {
  double I = 0.0;
  size_t N = X.size()-1;
  for ( size_t i = 0; i < N; i++ ) {
    I += V( X[i], t[i] ) * ( t[i+1] - t[i] );
  }
  return exp( -I ) * psi( X[N] );
}

int main( int argc, char* argv[] ) {
  size_t M = 100;
  size_t N = 10000;
  double t0, t1, x;
  Vector< double > t, G( M );
  Vector< Vector< double > > X( M );
 
  t0 = 0.0;
  t1 = 1.0;
  x = 1.0;
  
  for ( size_t i = 0; i < M; i++ ) {
    solveSTOCH( X[i], t, t0, t1, x, N, &u, &s );
    G.push_back( g( X[i], t ) );
  }
  
  cout << "E[g(X,t)|X_t = x] = " << sum( G ) / M << endl;
  
  return 0;
    
}