/*__________________________________________________________________________________________________
  autor: Pedro Guarderas
  email: ajusworkopensource@gmail.com
  date: 08-02-2015
  file: exo_4.cpp
  This program is free software; you can redistribute it and/or modify it under the 
  terms of the GNU General Public License as published by the Free Software Foundation; 
  either version 2 of the License, or (at your option) any later version.
  __________________________________________________________________________________________________
*/

#include <cmath>
#include <iostream>


using namespace std;


/*__________________________________________________________________________________________________
 * Cumulated normal distribution
 */
double pnorm( double x ) {
  double p = 0.2316419;
  double b1 = 0.319381530;
  double b2 = -0.356563782;
  double b3 = 1.781477937;
  double b4 = -1.821255978;
  double b5 = 1.330274429;
  double t = 1/(1 + p*x);
  
  double N = 1.0 - ( 1.0 / sqrt( 2.0 * 3.1415926535897932384626433 ) ) * exp( -x * x / 2.0 ) * 
  ( b1 * t + b2 * pow( t, 2 ) + b3 * pow( t, 3 ) + b4 * pow( t, 4 ) + b5 * pow( t, 5 ) );
  return N;
}

double pnorm( double x, double u, double s ) {
  return pnorm( ( x - u ) / s );
}

int main( int argc, char* argv[] ) {
  cout << pnorm( 0.0 ) << endl;
  cout << pnorm( 1.0, 1.0, 1.0 ) << endl;
  cout << pnorm( 3.0 ) << endl;
  cout << pnorm( 2.5, 1.0, 0.5 ) << endl;
  return 0;
}