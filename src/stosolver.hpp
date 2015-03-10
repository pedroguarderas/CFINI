/*__________________________________________________________________________________________________

  autor: Pedro Guarderas
  email: ajusworkopensource@gmail.com
  date: 09-03-2015
  file: stosolver.hpp
  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation;
  either version 2 of the License, or (at your option) any later version.
  __________________________________________________________________________________________________
*/

#ifndef STOSOLVER
#define STOSOLVER

#include <chrono>
#include <random>

#include "linear.hpp"


void solveSTOCH( Vector< double >& X, Vector< double >& t, 
                double t0, double t1, double x, size_t N,
                double (*u)( double, double ),
                double (*s)( double, double ) ) {
  typedef chrono::high_resolution_clock clock;
  
  size_t i;
  double dt;
  default_random_engine engine;
  normal_distribution< double > distribution( 0.0, 1.0 );
  clock::time_point beginning = clock::now();
  clock::duration duration;
  
  dt = ( t1 - t0 ) / ( N - 1.0 );
  X.push_back( x );
  t.push_back( t0 );
  
  for ( i = 1; i < N; i++ ) {
    t.push_back( t0 + i * dt );
    duration = clock::now() - beginning;
    engine.seed( duration.count() );
    X.push_back( X[i-1] + u( X[i-1], t[i-1] ) * dt + 
      sqrt( dt ) * s( X[i-1], t[i-1] ) * distribution( engine ) );
  }
}

#endif