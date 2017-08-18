#include <RcppArmadillo.h>

#include "TriDiagSolver.h"

using namespace Rcpp;

//[[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List DiffusionSolverES( const double& alpha,
                      const arma::colvec& I,
                      const arma::colvec& A,
                      const arma::colvec& B,
                      const arma::colvec& t,
                      const arma::colvec& x ) {
  
  int n, i;
  double dt, dxf, dxb, lambda, h;
  int Nt, Nx;
  
  Nt = t.size();
  Nx = x.size();
  arma::mat u( Nt, Nx );
  
  // Inicial condition
  for ( i = 0; i < Nx; i++ ) {
    u( 0, i ) = I( i );
  }
  
  // Boundary conditions
  for ( n = 0; n < Nt; n++ ) {
    u( n, 0 ) = A( n );
    u( n, Nt-1 ) = B( n );
  }
  
  // Solver
  // Euler explicit scheme
  for ( n = 0; n < Nt - 1; n++ ) {
    dt = t( n + 1 ) - t( n );
    for ( i = 1; i < Nx - 1; i++ ) {
      dxf = x( i + 1 ) - x( i );
      dxb = x( i ) - x( i - 1 );
      lambda = alpha * dt / ( dxf * dxb );
      u( n + 1, i ) = u( n, i ) + lambda * ( u( n, i + 1 ) - 2.0 * u( n, i ) + u( n, i - 1 ) );
    }
  }
  
  return List::create( Named( "u" ) = u,
                       Named( "t" ) = t,
                       Named( "x" ) = x );
}

// [[Rcpp::export]]
List DiffusionSolverCNS( const double& alpha,
                         const double& theta,
                         const arma::colvec& I,
                         const arma::colvec& A,
                         const arma::colvec& B,
                         const arma::colvec& t,
                         const arma::colvec& x ) {
  
  int n, i;
  double dt, dx, lambda, h;
  int Nt, Nx, nt, nx;
  Nt = t.size();
  Nx = x.size();
  
  nx = Nx - 1;
  nt = Nt - 1;
  
  arma::colvec a( Nx );
  arma::colvec b( Nx );
  arma::colvec d( Nx );
  arma::mat u( Nt, Nx );
  
  // Inicial condition
  for ( i = 0; i < Nx; i++ ) {
    u( 0, i ) = I( i );
  }
  
  // Boundary conditions
  for ( n = 0; n < Nt; n++ ) {
    u( n, 0 ) = A( n );
    u( n, Nx-1 ) = B( n );
  }
  
  // Solver
  // Crank-Nicolson scheme
  for ( n = 0; n < nt; n++ ) {
    dt = t( n + 1 ) - t( n );
    
    d( 0 ) = A( n );
    d( nx ) = B( n );
    
    dx = x( 1 ) - x( 0 );
    lambda = alpha * dt / ( dx * dx );
    a( 0 ) = -theta * lambda;
    b( 0 ) = 1.0 + 2.0 * theta * lambda;
    
    dx = x( nx ) - x( nx - 1 );
    lambda = alpha * dt / ( dx * dx );
    a( nx ) = -theta * lambda;
    b( nx ) = 1.0 + 2.0 * theta * lambda;
    
    for ( i = 1; i < nx; i++ ) {
      dx = x( i + 1 ) - x( i );
      lambda = alpha * dt / ( dx * dx );
      a( i ) = -theta * lambda;
      b( i ) = 1.0 + 2.0 * theta * lambda;
      
      d( i ) = u( n, i ) + ( 1.0 - theta ) * lambda * ( u( n, i + 1 ) - 2.0 * u( n, i ) + u( n, i - 1 ) );
    }
    
    TriDiagSolver( a, b, a, d );
    d( 0 ) = A( n );
    d( nx ) = B( n );
    
    for ( i = 0; i < Nx; i++ ) {
      u( n + 1, i ) = d( i );
    }

  }
  
  return List::create( Named( "u" ) = u,
                       Named( "t" ) = t,
                       Named( "x" ) = x );
}