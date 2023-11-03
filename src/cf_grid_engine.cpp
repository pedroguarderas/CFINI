
#include "cf_grid_engine.h"

//--------------------------------------------------------------------------------------------------
Eigen::VectorXd cf_uniform_grid( const double& a, const double& b, const double& n ) {
  double i;
  double h;
  Eigen::VectorXd X( static_cast< size_t >( n ) );
  
  h = ( b - a ) / ( n - 1.0 );
  
  for ( i = 0; i < n; i++ ) {
    X(i) = a + i * h;
  }
  return X;
}

//--------------------------------------------------------------------------------------------------
Eigen::VectorXd cf_adapt_grid( const double& l, const double& a, const double& b, const double& n,
                               const double& E ) {
  double i;
  Eigen::VectorXd X( static_cast< size_t >( n ) );
  double x, y, D, h;
  
  h = 1.0 / ( n - 1.0 );
  
  D = E * ( exp( l ) - 1.0 );
  x = ( a * E * exp( l ) - b * E ) / D;
  y = ( b - a ) / D;
  
  for ( i = 0; i < n; i++ ) {
    X(i) = x + y * exp( i * h * l );
  }
  return X;
}
