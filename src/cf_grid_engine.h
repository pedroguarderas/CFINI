#ifndef __CF_GRID_ENGINE__
#define __CF_GRID_ENGINE__

#include <RcppEigen.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

//--------------------------------------------------------------------------------------------------
//' @title Uniform grid
//' @description simple function for grid generation
//' @param a inferior value for the grid
//' @param b superior value for the grid
//' @param n number of points in the grid
//' @return A vector with the grid points.
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
Eigen::VectorXd cf_uniform_grid( const double& a, 
                                 const double& b, 
                                 const double& N );

//--------------------------------------------------------------------------------------------------
//' @title Exponential grid
//' @description Adaptive grid generator
//' @param l adaptation points
//' @param a inferior value for the grid
//' @param b superior value for the grid
//' @param n number of points in the grid
//' @param E refinement parameter
//' @return A vector with the adaptive grid points.
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
Eigen::VectorXd cf_adapt_grid( const double& l, 
                               const double& a,
                               const double& b,
                               const double& N,
                               const double& E );


#endif
