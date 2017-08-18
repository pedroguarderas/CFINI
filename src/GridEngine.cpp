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
    X(i) = a + i * h;
  }
  return X;
}

// [[Rcpp::export]]
arma::colvec GridExpAddapt( const double& l, const double& a, const double& b, const double& N ) {
  double i;
  arma::colvec X( N );
  double x, y, D, h;
  
  h = 1.0 / ( N - 1.0 );
  
  D = exp( l ) - 1.0;
  x = ( a * exp( l ) - b ) / D;
  y = ( b - a ) / D;
  
  for ( i = 0; i < N; i++ ) {
    X(i) = x + y * exp( i * h * l );
  }
  return X;
}
