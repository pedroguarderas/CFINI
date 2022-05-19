#ifndef __CF_BROWNIAN__
#define __CF_BROWNIAN__

#include <random>
#include <RcppEigen.h>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppEigen)]]

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
Eigen::MatrixXd cf_brownian( const int& d,
                             const Eigen::VectorXd& t );

//--------------------------------------------------------------------------------------------------
/*
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
Eigen::MatrixXd cf_stoch_solv( const int& d1,
                               const int& d2,
                               const Eigen::VectorXd& X0,
                               const Function& b,
                               const Function& s,
                               const Eigen::VectorXd& t );
*/
#endif
