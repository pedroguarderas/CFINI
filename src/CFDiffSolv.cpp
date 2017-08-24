#include <RcppArmadillo.h>

//[[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include "CFTriDiagSolv.h"

using namespace Rcpp;

//' @title Diffusion solver
//' @description Solver for standard diffusion problems
//' @param alpha
//' @param I
//' @param A
//' @param B
//' @param t
//' @param x
//' @return List with solution parameters
//' @note Diffusion solver for pricing options
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
List CFDiffSolvES( const double& alpha,
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

//' @title Diffusion solver
//' @description Solver for diffusion problems implemented with Crank-Nicolson scheme
//' @param alpha
//' @param theta
//' @param I
//' @param A
//' @param B
//' @param t
//' @param x
//' @return List with solution parameters
//' @note Diffusion solver for pricing options
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
List CFDiffSolvCNS( const double& alpha,
                    const double& theta,
                    const arma::colvec& I,
                    const arma::colvec& A,
                    const arma::colvec& B,
                    const arma::colvec& t,
                    const arma::colvec& x ) {
  
  int n, i;
  double dt, dxf, dxb, lambdaf, lambdab, h;
  int Nt, Nx, nt, nx;
  Nt = t.size();
  Nx = x.size();
  
  nx = Nx - 1;
  nt = Nt - 1;
  
  arma::colvec a( Nx );
  arma::colvec b( Nx );
  arma::colvec c( Nx );
  arma::colvec d( Nx );
  arma::mat u( Nt, Nx );
  
  // Inicial condition
  for ( i = 0; i < Nx; i++ ) {
    u( 0, i ) = I( i );
  }
  
  // Boundary conditions
  for ( n = 0; n < Nt; n++ ) {
    u( n, 0 ) = A( n );
    u( n, nx ) = B( n );
  }
  
  // Solver
  // Crank-Nicolson scheme
  for ( n = 0; n < nt; n++ ) {
    dt = t( n + 1 ) - t( n );
    
    d( 0 ) = A( n );
    d( nx ) = B( n );
    
    dxf = x( 1 ) - x( 0 );
    lambdaf = alpha * dt / ( dxf * dxf );

    dxb = x( nx ) - x( nx - 1 );
    lambdab = alpha * dt / ( dxb * dxf );
    
    a( 0 ) = -theta * lambdab;
    a( nx ) = -theta * lambdab;
    b( 0 ) = 1.0 + theta * ( lambdaf + lambdab );
    b( nx ) = 1.0 + theta * ( lambdaf + lambdab );
    c( 0 ) = -theta * lambdaf;
    c( nx ) = -theta * lambdaf;
    
    for ( i = 1; i < nx; i++ ) {
      dxf = x( i + 1 ) - x( i );
      dxb = x( i ) - x( i - 1 );
      lambdab = alpha * dt / ( dxb * dxf );
      lambdaf = alpha * dt / ( dxf * dxf );
      
      a( i ) = -theta * lambdab;
      b( i ) = 1.0 + theta * ( lambdaf + lambdab );
      c( i ) = -theta * lambdaf;
      
      d( i ) = u( n, i ) + ( 1.0 - theta ) * 
        ( lambdaf * ( u( n, i + 1 ) - u( n, i ) ) - lambdab * ( u( n, i ) - u( n, i - 1 ) ) ) ;
    }
    
    CFTriDiagSolv( a, b, c, d );
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

//' @title Black-Scholes solver
//' @description Solver for Black-Scholes models implemented with Crank-Nicolson scheme
//' @param sigma
//' @param rate
//' @param theta
//' @param I
//' @param A
//' @param B
//' @param t
//' @param x
//' @return List with solution parameters
//' @note pricing options
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
List CFBlackScholesSolvCNS( const double& sigma,
                            const double& rate,
                            const double& theta,
                            const arma::colvec& I,
                            const arma::colvec& A,
                            const arma::colvec& B,
                            const arma::colvec& t,
                            const arma::colvec& x ) {
  
  int n, i;
  double dt, dxf, dxb, lambdaf, lambdab, h;
  double alpha, beta;
  int Nt, Nx, nt, nx;
  Nt = t.size();
  Nx = x.size();
  
  alpha = ( sigma * sigma - 2.0 * rate ) / ( 2.0 * sigma * sigma );
  beta = ( 2.0 * rate + sigma * sigma ) / ( 2.0 * sigma * sigma );
  beta = -beta * beta;
  
  nx = Nx - 1;
  nt = Nt - 1;
  
  arma::colvec a( Nx );
  arma::colvec b( Nx );
  arma::colvec c( Nx );
  arma::colvec d( Nx );
  
  arma::mat u( Nt, Nx );
  arma::mat v( Nt, Nx );
  
  // Inicial condition
  for ( i = 0; i < Nx; i++ ) {
    u( 0, i ) = I( i );
  }
  
  // Boundary conditions
  for ( n = 0; n < Nt; n++ ) {
    u( n, 0 ) = A( n );
    u( n, nx ) = B( n );
  }
  
  // Solver
  // Crank-Nicolson scheme
  for ( n = 0; n < nt; n++ ) {
    dt = t( n + 1 ) - t( n );
    
    d( 0 ) = A( n );
    d( nx ) = B( n );
    
    dxf = x( 1 ) - x( 0 );
    lambdaf = sigma * dt / ( dxf * dxf );
    
    dxb = x( nx ) - x( nx - 1 );
    lambdab = sigma * dt / ( dxb * dxf );
    
    a( 0 ) = -theta * lambdab;
    a( nx ) = -theta * lambdab;
    b( 0 ) = 1.0 + theta * ( lambdaf + lambdab );
    b( nx ) = 1.0 + theta * ( lambdaf + lambdab );
    c( 0 ) = -theta * lambdaf;
    c( nx ) = -theta * lambdaf;
    
    for ( i = 1; i < nx; i++ ) {
      dxf = x( i + 1 ) - x( i );
      dxb = x( i ) - x( i - 1 );
      lambdab = sigma * dt / ( dxb * dxf );
      lambdaf = sigma * dt / ( dxf * dxf );
      
      a( i ) = -theta * lambdab;
      b( i ) = 1.0 + theta * ( lambdaf + lambdab );
      c( i ) = -theta * lambdaf;
      
      d( i ) = u( n, i ) + ( 1.0 - theta ) * 
        ( lambdaf * ( u( n, i + 1 ) - u( n, i ) ) - lambdab * ( u( n, i ) - u( n, i - 1 ) ) ) ;
    }
    
    CFTriDiagSolv( a, b, c, d );
    d( 0 ) = A( n );
    d( nx ) = B( n );
    
    for ( i = 0; i < Nx; i++ ) {
      u( n + 1, i ) = d( i );
    }
    
  }
  
  for ( i = 0; i < Nx; i++ ) {
    for ( n = 0; n < Nt; n++ ) {
      v( n, i ) = exp( alpha * x( i ) ) * exp( beta * ( t( nt ) - t( nt - n ) ) ) * u( nt - n, i );
    }
    d( i ) = exp( x(i) );
  }
  
  return List::create( Named( "u" ) = v,
                       Named( "t" ) = t,
                       Named( "x" ) = d );
}
