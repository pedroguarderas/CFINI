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
#include "stosolver.hpp"


using namespace std;


double u( double x, double t ) {
  return sin( 3.14159 * x );
}

double s( double x, double t ) {
  return 0.5*sqrt(x);
}

int main( int argc, char* argv[] ) {
  size_t N;
  double t0, t1, x;
  Vector< double > t, X;
 
  N = 10000;
  t0 = 0.0;
  t1 = 1.0;
  x = 1.0;
  
  solveSTOCH( X, t, t0, t1, x, N, &u, &s );

  ofstream file;
  file.open ("ito_process.txt");
  for ( size_t i = 0; i < N; i++ ) {
    file << t[i] << "\t" << X[i] << endl;
  }
  file.close();
  
  return 0;
    
}