#ifndef __CF_THIELE_SOLV__
#define __CF_THIELE_SOLV__

#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//--------------------------------------------------------------------------------------------------
//' @title Thiele equations solver
//' @description Solver for Thiele equation and computation of Mathematical reserves. The solver 
//' is implemented with a Crank-Nicolson algorithm.
//' @param t time grid, could be adapted
//' @param V0 initial mathematical reserve
//' @param b fixed benefits
//' @param B transition benefits
//' @return Return a list with the mathematical reserve, the time grid.
//' @note The solver is implemented to compute the solution forwards in time.
//' @author Pedro Guarderas
//' @useDynLib CFINI
//' @importFrom Rcpp sourceCpp
//' @exportPattern("^[[:alpha:]]+")
//' @export
// [[Rcpp::export]]
List cf_thiele_solv( const arma::colvec& t,
                     const arma::colvec& V0,
                     const arma::mat& b,
                     const arma::cube& B,
                     const double& theta = 0.5 );

#endif
