
#include <algorithm>
#include "cf_tri_diag_solv.h"

//--------------------------------------------------------------------------------------------------
void cf_tri_diag_solv( Eigen::VectorXd& a,
                       Eigen::VectorXd& b,
                       Eigen::VectorXd& c,
                       Eigen::VectorXd& d ) {
  
  int N = d.size();
  int i;
  
  N--; 
  c(0) /= b(0);
  d(0) /= b(0);
  
  // Forward recursion
  for ( i = 1; i < N; i++ ) {
    c(i) /= b( i ) - a( i ) * c( i - 1 );
    d( i ) = ( d( i ) - a( i ) * d( i - 1 ) ) / ( b( i ) - a( i ) * c( i - 1 ) );
  }
  
  d( N ) = ( d( N ) - a( N ) * d( N - 1 ) ) / ( b( N ) - a( N ) * c( N - 1 ) );
  
  // Backsubstitution
  for ( i = N; i-- > 0; ) {
    d( i ) -= c( i ) * d( i + 1 );
  }
  
}

//--------------------------------------------------------------------------------------------------
List cf_sor_solv( const Eigen::VectorXd& u0,
                  const Eigen::MatrixXd& A, 
                  const Eigen::VectorXd& b,
                  const Eigen::VectorXd& c,
                  const double& w,
                  const int& n,
                  const double& e ) {
  
  int i, j, k, ek;
  int d = b.size();
  double v;
  
  Eigen::VectorXd u( d ), u1( d );
  u = u0;
  
  while( k <= n || ek > e ) {
    u1 = u;
    for ( i = 0; i < d; i++ ) {
      v = b( i ) + ( 1 / w  - 1 ) * A( i, i ) * u( i );
      for ( j > i; i < d; i++ ) {
        v = v - A( i, j ) * u( j ) - A( j, i ) * u( j );
      }
      u( i ) = ( w / A( i, i ) ) * v;
      u( i ) = std::max( u( i ), c( i ) );  
    }
    ek = ( u - u1 ).norm() / u.norm();
    k++;
  }
  
  return List::create( Named( "u" ) = u, Named( "k" ) = k , Named( "e" ) = e );
  
}

