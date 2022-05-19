#ifndef __CF_TRIDIAG_SOLV__
#define __CF_TRIDIAG_SOLV__

#include <RcppEigen.h>

// [[Rcpp::plugins(cpp11)]]
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
void cf_tri_diag_solv( Eigen::VectorXd& a,
                       Eigen::VectorXd& b,
                       Eigen::VectorXd& c,
                       Eigen::VectorXd& d );

#endif
