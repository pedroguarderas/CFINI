/*__________________________________________________________________________________________________
  autor: Pedro Guarderas
  email: ajusworkopensource@gmail.com
  date: 20-07-2015
  file: exo_12.cpp
  This program is free software; you can redistribute it and/or modify it under the 
  terms of the GNU General Public License as published by the Free Software Foundation; 
  either version 2 of the License, or (at your option) any later version.
  __________________________________________________________________________________________________
*/

#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include "linear.hpp"

using namespace std;


double Ue( double x, double y ) {
    return( ( 1.0 - x*x - y*y ) / 4.0 );
}


/*__________________________________________________________________________________________________
 * Normal random number generation
 */
int main( int argc, char* argv[] ) {
  default_random_engine engine;
  normal_distribution< double > distribution( 0.0, 1.0 );
  Vector< double > X, Y;
  Vector< double > x(2), B(2), W(2);
  size_t n = 20;
  Matrix< double > u( n, n );
  Matrix< double > ue( n, n );
  double dt = 0.001, d;
  double E, M = 100, tau;
  size_t i, j, k;
  
  d = ( 1.0 + 1.0 ) / ( n - 1 );
  for ( size_t i = 0; i < n; i++ ) {
    X.push_back( -1.0 + d * i );
    Y.push_back( -1.0 + d * i );
  }
  // shared(dt,u,ue,M,n,X,Y)
//   #pragma omp parallel for shared(X,Y,u,ue,n,M,dt) private(k,B,E,tau) collapse(2)
  for (  i = 0; i < n; i++ ) {
      for (  j = 0; j < n; j++ ) {

          x[0] = X[i]; x[1] = Y[j];

          if ( norm( x ) <  1.0 ) {
              
            E = 0.0;
            for ( k = 0; k < M; k++ ) {
              B = x;
              tau = 0;
              while ( norm( B ) < 1 ) {
                  W[0] = distribution( engine );
                  W[1] = distribution( engine );
                  B = B + sqrt( dt ) * W;
                  tau++;
              }
             
              E = E + dt * tau;
            }
            
            u(i,j) = 0.5 * E / M;
            ue(i,j) = Ue( x[0], x[1] );
          } else {
            u(i,j) = 0.0;
            ue(i,j) = 0.0;
          }
      }
  }
  
  ofstream file;
  file.open ("heat_stochastic_solution.txt");
  file << u << endl;
  file.close();
  
  file.open ("heat_stochastic_solution_exact.txt");
  file << ue << endl;
  file.close();
  
  return 0;
}