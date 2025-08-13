
#include "cf_grid_engine.h"

//--------------------------------------------------------------------------------------------------
Eigen::VectorXd cf_uniform_grid( const double& a, const double& b, const unsigned int& n ) {
  
  unsigned int i;
  double h;
  Eigen::VectorXd X( n );
  
  if ( n == 1 ) {
    
    X( 0 ) = a;
    
  } else if ( n > 1 ) {
    
    h = ( b - a ) / ( n - 1.0 );
    
    #pragma omp parallel for
    for ( i = 0; i < n; i++ ) {
      
      X( i ) = a + i * h;
      
    }
    
  }
  
  return X;
}

//--------------------------------------------------------------------------------------------------
Eigen::VectorXd cf_adapt_grid( const double& l, const double& a, const double& b, 
                               const unsigned int& n, const double& E ) {
  unsigned int i;
  Eigen::VectorXd X( n );
  double x, y, D, h;
  
  h = 1.0 / ( n - 1.0 );
  
  D = E * ( exp( l ) - 1.0 );
  x = ( a * E * exp( l ) - b * E ) / D;
  y = ( b - a ) / D;
  
  #pragma omp parallel for
  for ( i = 0; i < n; i++ ) {
    
    X( i ) = x + y * exp( i * h * l );
    
  }
  
  return X;
}
