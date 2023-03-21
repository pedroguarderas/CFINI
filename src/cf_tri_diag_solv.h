#ifndef __CF_TRI_DIAG_SOLV__
#define __CF_TRI_DIAG_SOLV__

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



//--------------------------------------------------------------------------------------------------
//' @title PSOR algorithm
//' @description Projected successive over-relaxation
//' @param u0 initial guest of the solution
//' @param A matrix determining the variational inequality
//' @param b
//' @param c lower constraint for the solution
//' @param w weight
//' @param n maximal number of iterations
//' @param e relative error for the solution improvement
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
List cf_psor_solv( const Eigen::VectorXd& u0,
                   const Eigen::MatrixXd& A, 
                   const Eigen::VectorXd& b,
                   const Eigen::VectorXd& c,
                   const double& w,
                   const int& n,
                   const double& e );

#endif

