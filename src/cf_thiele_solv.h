#ifndef __CF_THIELE_SOLV__
#define __CF_THIELE_SOLV__

#include <RcppEigen.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

//--------------------------------------------------------------------------------------------------
/*
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
List cf_thiele_solv( const Eigen::VectorXd& t,
                     const Eigen::VectorXd& V0,
                     const Eigen::MatrixXd& b,
                     const Eigen::MatrixXd& B,
                     const double& theta = 0.5 );
*/
#endif
