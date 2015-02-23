/*__________________________________________________________________________________________________

  autor: Andrés Lopez
  email:
  date: 23-02-2015
  file: exo_5.cpp
  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation;
  either version 2 of the License, or (at your option) any later version.
  __________________________________________________________________________________________________
*/


#include <iostream>
#include <vector>

using namespace std;


void solve( const vector< double > a, const vector< double > b,
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

int main() {
	vector< double > a = { 0, -1, -1, -1 };
	vector< double > b = { 4,  4,  4,  4 };
	vector< double > c = { -1, -1, -1,  0 };
	vector< double > d = { 5,  5, 10, 23 };
	
  solve( a, b, c, d );

	for( size_t i = 0; i < a.size(); i++ ) {
		cout << d[i] << endl;
	}

	return 0;
}




