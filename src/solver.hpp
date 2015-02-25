/*__________________________________________________________________________________________________

  autor: Andrés López
  email: andrelo1d@hotmail.com
  date: 23-02-2015
  file: exo_5.cpp
  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation;
  either version 2 of the License, or (at your option) any later version.
  __________________________________________________________________________________________________
*/

#include <vector>

using namespace std;

void solveTDS( const vector< double > a, const vector< double > b,
            vector< double >& c, vector< double >& d ) {

	size_t n = a.size();
	c[0] /= b[0];
	d[0] /= b[0];
	n--;

	for( size_t i = 1; i < n; i++ ) {
		c[i] /= b[i] - a[i] * c[i - 1];
		d[i] = ( d[i] - a[i] * d[i - 1] ) / ( b[i] - a[i] * c[i - 1] );
	}

	d[n] = ( d[n] - a[n] * d[n - 1] ) / ( b[n] - a[n] * c[n - 1] );

	for( size_t i = n; i-- > 0; ) {
		d[i] -= c[i] * d[i + 1];
	}
}

void solveTDSnd( const vector< double > a, const vector< double > b,
            vector< double >& C, vector< double >& d ) {

  vector< double > c = C;
  size_t n = a.size();
  c[0] /= b[0];
  d[0] /= b[0];
  n--;

  for( size_t i = 1; i < n; i++ ) {
    c[i] /= b[i] - a[i] * c[i - 1];
    d[i] = ( d[i] - a[i] * d[i - 1] ) / ( b[i] - a[i] * c[i - 1] );
  }

  d[n] = ( d[n] - a[n] * d[n - 1] ) / ( b[n] - a[n] * c[n - 1] );

  for( size_t i = n; i-- > 0; ) {
    d[i] -= c[i] * d[i + 1];
  }
}

void multTDS( const vector< double >& a, const vector< double >& b,
              const vector< double >& c, vector< double >& d ) {
  size_t n = b.size()-1;
  d[0] = b[0] * d[0] + c[0] * d[1];
  for ( size_t i = 1; i < n; i++ ) {
    d[i] = a[i-1] * d[i-1] + b[i] * d[i] + c[i+1] * d[i+1];
  }
  d[n] = a[n] * d[n-1] + c[n] * d[n];
}
