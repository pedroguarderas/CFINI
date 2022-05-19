#ifndef __CF_TRIDIAG_SOLV__
#define __CF_TRIDIAG_SOLV__

#include <RcppEigen.h>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

//--------------------------------------------------------------------------------------------------
//' @title Tridiagonal solver
//' @description Solver tridiagonal matrices.
//' @param alpha
//' @param a lower diagonal
//' @param b diagonal
//' @param c upper diagonal
//' @param d image of the solution vector, which over written with the solution
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
void cf_tri_diag_solv( arma::colvec& a,
                       arma::colvec& b,
                       arma::colvec& c,
                       arma::colvec& d );

#endif
