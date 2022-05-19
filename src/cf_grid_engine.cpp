
#include "cf_grid_engine.h"

//--------------------------------------------------------------------------------------------------
Eigen::VectorXd cf_uniform_grid( const double& a, const double& b, const double& N ) {
  double i;
  double h;
  Eigen::VectorXd X( static_cast< size_t >( N ) );
  
  h = ( b - a ) / ( N - 1.0 );
  
  for ( i = 0; i < N; i++ ) {
    X(i) = a + i * h;
  }
  return X;
}

//--------------------------------------------------------------------------------------------------
Eigen::VectorXd cf_adapt_grid( const double& l, const double& a, const double& b, const double& N,
                               const double& E ) {
  double i;
  Eigen::VectorXd X( static_cast< size_t >( N ) );
  double x, y, D, h;
  
  h = 1.0 / ( N - 1.0 );
  
  D = E * ( exp( l ) - 1.0 );
  x = ( a * E * exp( l ) - b * E ) / D;
  y = ( b - a ) / D;
  
  for ( i = 0; i < N; i++ ) {
    X(i) = x + y * exp( i * h * l );
  }
  return X;
}
