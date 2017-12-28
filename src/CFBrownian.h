#ifndef __CF_BROWNIAN__
#define __CF_BROWNIAN__

#include <random>
#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//--------------------------------------------------------------------------------------------------
//' @title Brownian motion
//' @description Simulate d-dimensional Browninan motion
//' @param d Dimension
//' @param t Time grid
//' @return List with solution parameters
//' @note Diffusion solver for pricing options
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
arma::mat CFRBrownian( const int& d,
                       const arma::colvec& t );

//--------------------------------------------------------------------------------------------------
//' @title Brownian motion
//' @description Simulate d-dimensional Brownian motion
//' @param d1
//' @param d2
//' @param X0
//' @param b
//' @param s
//' @param t
//' @return Solution to the stochastic differential equation
//' @note Diffusion solver for pricing options
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
arma::mat CFStochSolv( const int& d1,
                       const int& d2,
                       const arma::colvec& X0,
                       const Function& b,
                       const Function& s,
                       const arma::colvec& t );

#endif


