/*__________________________________________________________________________________________________
  autor: Pedro Guarderas
  email: ajusworkopensource@gmail.com
  date: 28-01-2014
  file: exo_3.cpp
  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation;
  either version 2 of the License, or (at your option) any later version.
  __________________________________________________________________________________________________
*/

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

using namespace std;


/*__________________________________________________________________________________________________
 * Brownian motion simulation
 */
double CallPayoff( double S, double K ) {
  double X = 0.0;
  if ( S > K ) {
     X = S - K;
  }
  return X;
}

double S( double S0, double u, double d, int n, int i ) {
  return S0 * pow( 1 + u, i ) * pow( 1 + d, n - i );
}

double RiskNeutralProb( double u, double d, double r ) {
  return ( r - d ) / ( u - d );
}

double PriceCRR( double S0, double u, double d, double r, int T, double K ) {
    double p;
    vector< double > C;
    
    p = RiskNeutralProb( u, d, r );
    for ( int i = 0; i <= T; i++ ) {
      C.push_back( CallPayoff( S( S0, u, d, T, i ), K ) );
    }
    
    for ( int n = T - 1; n >= 0; n-- ) {
      for ( int i = 0; i <= n; i++ ) {
        C[i] = ( p * C[i+1] + ( 1 - p ) * C[i] ) / ( 1 + r );
      }
    }
    return C[0];
}


int main( int argc, char* argv[] ) {
  double S0, u, d, r, K;
  int T;
  S0 = 8.0;
  u = 0.05;
  d = -0.06;
  r = 0.01;
  K = 8.1;
  T = 1;
  
  // writing in external file
  cout << "Risk Neutral Probability = " << RiskNeutralProb( u, d, r ) << endl;
  cout << "European call option price = " << PriceCRR( S0, u, d, r, T, K ) << endl;
 
  return 0;
}