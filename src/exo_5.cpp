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


#include <iostream>
#include <vector>

#include "solver.hpp"

using namespace std;


int main() {
	vector< double > a = { -1, -1, -1 };
	vector< double > b = { 4,  4,  4,  4 };
	vector< double > c = { -1, -1, -1 };
	vector< double > d = { 5,  5, 10, 23 };
	
	solveTDS( a, b, c, d );

	for( size_t i = 0; i < d.size(); i++ ) {
		cout << d[i] << endl;
	}

	return 0;
}




