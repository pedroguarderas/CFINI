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

/*__________________________________________________________________________________________________
 * Solving second order differential equation
 */
double f( double t ) {
	return( 100 * sin( 3.14159 * t ) * cos( 20 * 3.14159 * t ) + 0.5 * t );
}

int main() {
	double n = 10000;
	double T0 = 0.0, T1 = 2.0, h;
	double alpha, beta;
	vector< double > t( n + 2 ), F( n ), a( n ), b( n ), c( n );

	alpha = 0.5;
	beta = 1.0;
	h = ( T1 - T0 ) / ( n + 1 );

	t[0] = T0;
	for( size_t i = 0; i < n; i++ ) {
		t[i + 1] = T0 + ( i + 1 ) * h;
		F[i] = f( t[i + 1] );
		a[i] = 1.0 / pow( h, 2 );
		b[i] = -2.0 / pow( h, 2 );
		c[i] = 1.0 / pow( h, 2 );
	}
	t[n + 1] = T1;
	F[0] -= alpha / pow( h, 2 );
	F[n - 1] -= beta / pow( h, 2 );

	solveTDS( a, b, c, F );

	ofstream file;
	file.open( "output.txt" );
	file << t[0] << "\t" << alpha << endl;
	for( size_t i = 0; i < n; i++ ) {
		file << t[i + 1] << "\t" << F[i] << endl;
	}
	file << t[n + 1] << "\t" << beta << endl;
	file.close();

	return 0;
}

/* Gráfico en R
 s<-read.table( '[path]/CFINI/src/output.txt')
 plot( s[,1], s[,2], col = 'red', cex = 0.5, type = 'l', xlab = 't', ylab = 'u' )
 */