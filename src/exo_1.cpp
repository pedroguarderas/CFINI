/*__________________________________________________________________________________________________
  autor: Pedro Guarderas
  email: ajusworkopensource@gmail.com
  date: 28-01-2015
  file: exo_1.cpp
  This program is free software; you can redistribute it and/or modify it under the 
  terms of the GNU General Public License as published by the Free Software Foundation; 
  either version 2 of the License, or (at your option) any later version.
  __________________________________________________________________________________________________
*/

#include <fstream>
#include <iostream>
#include <random>
#include <vector>

using namespace std;


/*__________________________________________________________________________________________________
 * Normal random number generation
 */
int main( int argc, char* argv[] ) {
  default_random_engine engine;
  normal_distribution< double > distribution( 0.0, 1.0 );
  vector< double > v;
  size_t n = 100;
  
  for ( size_t i = 0; i < n; i++ ) {
    v.push_back( distribution( engine ) );
  }
  
  for ( size_t i = 0; i < n; i++ ) {
    cout << v[i] << " ";
  }
  cout << endl;
  
  return 0;
}
