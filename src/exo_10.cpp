/*__________________________________________________________________________________________________

  autor: Pedro Guarderas
  email: ajusworkopensource@gmail.com
  date: 09-03-2015
  file: exo_10.cpp
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

using namespace std;

double u( double x, double t ) {
  return sin( 3.14159 * x );
}

double s( double x, double t ) {
  return 0.25;
}

int main( int argc, char* argv[] ) {
  size_t N;
  double t0, t1, dt, x;
  Vector< double > t, X;
  default_random_engine engine;
  normal_distribution< double > distribution( 0.0, 1.0 );
  
  typedef chrono::high_resolution_clock clock;
  clock::time_point beginning = clock::now();
  clock::duration duration;
  
  N = 10000;
  t0 = 0.0;
  t1 = 1.0;
  x = 1.0;
  
  dt = ( t1 - t0 ) / ( N - 1.0 );
  X.push_back( x );
  t.push_back( t0 );
  
  for ( size_t i = 1; i < N; i++ ) {
    t.push_back( t0 + i * dt );
    duration = clock::now() - beginning;
    engine.seed( duration.count() );
    X.push_back( X[i-1] + u( X[i-1], t[i-1] ) * dt + 
      sqrt( dt ) * s( X[i-1], t[i-1] ) * distribution( engine ) );
  }

  ofstream file;
  file.open ("ito_process.txt");
  for ( size_t i = 0; i < N; i++ ) {
    file << t[i] << "\t" << X[i] << endl;
  }
  file.close();
  
  return 0;
    
}