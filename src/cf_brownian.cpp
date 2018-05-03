
#include "cf_brownian.h"

//--------------------------------------------------------------------------------------------------
arma::mat cf_brownian( const int& d,
                       const arma::colvec& t ) {
  
  int i, j, n;
  double sdt;
  n = t.size();
  
  arma::mat W( n, d );
  std::random_device engine;
  std::normal_distribution< double > Normal( 0, 1.0 );
  
  for ( j = 0; j < d; j++ ) {
    W( 0, j ) = 0.0;
  }
  
  for ( i = 1; i < n; i++ ) {
    sdt = sqrt( t( i ) - t( i - 1 ) );
    for ( j = 0; j < d; j++ ) {
      W( i, j ) = W( i - 1, j ) + sdt * Normal( engine );
    }
  }
  
  return( W );
}

//--------------------------------------------------------------------------------------------------
arma::mat cf_stoch_solv( const int& d1,
                         const int& d2,
                         const arma::colvec& X0,
                         const Function& b,
                         const Function& s,
                         const arma::colvec& t ) {
  int i, k, n;
  double dt, sdt;
  std::random_device engine;
  std::normal_distribution< double > Normal( 0, 1.0 );
  
  n = t.size();
  
  arma::mat X( n, d1 );
  arma::mat S( d1, d2 );
  arma::colvec B( d1 ), W( d2 ), dW( d2 );
  
  X.row( 0 ) = X0.t();
  W.zeros();
  
  for ( i = 1; i < n; i++ ) {
    dt = t( i ) - t( i - 1 );
    sdt = sqrt( dt );
    
    B = as< arma::colvec >( b( t( i ), X.row( i - 1 ).t() ) );
    S = as< arma::mat >( s( t( i ), X.row( i - 1 ).t() ) );
    
    for ( k = 0; k < d2; k++ ) {
      dW( k ) = W( k ) + sdt * Normal( engine );
    }
    B = dt * B  + S * dW;
    X.row( i ) = X.row( i - 1 ) +  B.t();
    W = dW;
  }
  
  return( X );
}

