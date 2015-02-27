/*__________________________________________________________________________________________________

  autor: Pedro Guarderas
  email: ajusworkopensource@gmail.com
  date: 26-02-2015
  file: exo_8.cpp
  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation;
  either version 2 of the License, or (at your option) any later version.
  __________________________________________________________________________________________________
*/

#include <fstream>
#include <vector>
#include <cmath>

#include "solver.hpp"

using namespace std;

double f( double x ) {
  return( max( x - 20.0, 0.0 ) );
//   return( sin( sx ) );
}

double bd( double t ) {
  return( 80.0 );
}

double sigma( double t, double x ) {
  return( 1.0 );
}

double mu( double t, double x ) {
  return( 2.0 );
}

int main() {
  double n = 600;
  double m = 600;
  double T = 2.0, X = 100.0;
  double h, k;
  double r;
  
  vector< double > t( n ), a( n ), b( n ), c( n ); 
  vector< double > x( m ), s( m );
  vector< vector< double > > S;
  
  h = T / ( n - 1.0 );
  k = X / ( m - 1.0 );
  
  r = 1.05;
  
  for ( size_t j = 0; j < m; j++ ) {
    x[j] = j * k;
    s[j] = f( x[j] ); // Terminal condition
  }
 
  for ( size_t i = 0; i < n; i++ ) {
    t[i] = T - i * h; // Solving backward in time
    S.push_back( s );
    for ( size_t j = 0; j < m; j++ ) {
      a[j] = 0.5 * h * ( mu( t[i], x[j] ) / k - pow( sigma( t[i], x[j] ) / k, 2.0 ) ); 
      b[j] = 1 + h * ( r + pow( sigma( t[i], x[j] ) / k, 2.0 ) );
      c[j] = -0.5 * h * ( mu( t[i], x[j] ) / k + pow( sigma( t[i], x[j] ) / k, 2.0 ) );
    }
    solveTDS( a, b, c, s );
    s[m-1] = bd( t[i] ); // Including boundary condition
  } 

  ofstream file;
  file.open ("solution.txt");
  for ( size_t i = 0; i < n; i++ ) {
    for ( size_t j = 0; j < m; j++ ) {
      file << S[i][j] << "\t";
    }
    file << endl;
  }
  file.close();

  return 0;
}