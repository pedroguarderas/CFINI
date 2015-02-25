/*__________________________________________________________________________________________________

  autor: Pedro Guarderas
  email: ajusworkopensource@gmail.com
  date: 25-02-2015
  file: exo_6.cpp
  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation;
  either version 2 of the License, or (at your option) any later version.
  __________________________________________________________________________________________________
*/

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

#include "solver.hpp"

using namespace std;

double f( double t ) {
   return( sin(3.14159 * t ) * cos( 20 * 3.14159 * t )  ); 
}

int main() {
	double n = 1000;
	double T = 2.0, h;
	vector< double > t( n ), F( n ), a( n ), b( n ), c( n ); 
	h = T / ( n - 1 );
	
	for ( size_t i = 0; i < n; i++ ) {
	  t[i] = i * h;
	  F[i] = f( t[i] );
	  a[i] = 1.0 / pow( h, 2 );
	  b[i] = -2.0 / pow( h, 2 );
	  c[i] = 1.0 / pow( h, 2 );
	}

	solveTDS( a, b, c, F );

	for( size_t i = 0; i < n; i++ ) {
		cout << t[i] << "\t" << F[i] << endl;
	}

	return 0;
}

/* Gráfico en R
 * s<-read.table( '/home/drew/Development/CFINI/src/output.txt')
   plot( s[,1], s[,2], col = 'red', cex = 0.5, type = 'l' ) */