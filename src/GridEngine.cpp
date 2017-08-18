#include <RcppArmadillo.h>

using namespace Rcpp;

//[[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::colvec GridUniform( const double& a, const double& b, const double& N ) {
  double i;
  arma::colvec X( N );
  double h;
  
  h = ( b - a ) / ( N - 1.0 );
  
  for ( i = 0; i < N; i++ ) {
    X(i) = i * h;
  }
  return X;
}
