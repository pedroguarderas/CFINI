#include <RcppArmadillo.h>

using namespace Rcpp;

//[[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double GridEngine( const double a, const double b, const double N ) {
  double i;
  // arma::colvec X( N );
  double h;
  
  h = ( b - a ) / ( N - 1.0 );
  
  // for ( i = 0; i < N; i++ ) {
  //   X(i) = i * h;
  // }
  return h;
}

// [[Rcpp::export]]
List DiffusionSolver( const double& alpha,
                      const arma::colvec& I,
                      const arma::colvec& A,
                      const arma::colvec& B,
                      const double& t0, const double& t1, const double& Nt,
                      const double& x0, const double& x1, const double& Nx ) {

  int n, i;
  double dt, dx, lambda, h;

  arma::colvec x( Nx );
  arma::colvec t( Nt );
  arma::mat u( Nt, Nx );

  x = GridEngine( x0, x1, Nx );
  t = GridEngine( t0, t1, Nt );

  // Inicial condition
  for ( i = 0; i < Nx; i++ ) {
    u( 0, i ) = I( i );
  }

  // Boundary conditions
  for ( n = 0; n < Nt; n++ ) {
    u( n, 0 ) = A( n );
    u( n, Nt-1 ) = B( n );
  }

  // Explicit solver
  for ( n = 0; n < Nt - 1; n++ ) {
    dt = t( n + 1 ) - t( n );
    for ( i = 1; i < Nx - 1; i++ ) {
      dx = x( i + 1 ) - x( i );
      lambda = alpha * dt / ( dx * dx );
      u( n + 1, i ) = u( n, i ) + lambda * ( u( n, i + 1 ) - 2 * u( n, i ) + u( n, i - 1 ) );
    }
  }

  return List::create( Named( "u" ) = u,
                       Named( "t" ) = t,
                       Named( "x" ) = x );
}