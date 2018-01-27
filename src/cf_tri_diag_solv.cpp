
#include "cf_tri_diag_solv.h"

//--------------------------------------------------------------------------------------------------
void CFTriDiagSolv( arma::colvec& a,
                    arma::colvec& b,
                    arma::colvec& c,
                    arma::colvec& d ) {
  
  double mu;
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