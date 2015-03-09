/*__________________________________________________________________________________________________

  autor: Pedro Guarderas
  email: ajusworkopensource@gmail.com
  date: 09-03-2015
  file: exo_9.cpp
  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation;
  either version 2 of the License, or (at your option) any later version.
  __________________________________________________________________________________________________
*/

#include <iostream>
#include <random>
#include "linear.hpp"

using namespace std;


int main( int argc, char* argv[] ) {
  size_t n = 10;
  Vector< double > v, w, x;
  default_random_engine engine;
  normal_distribution< double > distribution( 0.0, 1.0 );
  
  for ( size_t i = 0; i < n; i++ ) {
    v.push_back( distribution( engine ) );
    w.push_back( distribution( engine ) );
  }
  
  x = v + w;
  cout << v << endl;
  cout << w << endl;
  cout << x << endl;
  return 0;
}