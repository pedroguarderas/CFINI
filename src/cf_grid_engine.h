#ifndef __CF_GRID_ENGINE__
#define __CF_GRID_ENGINE__

#include <RcppArmadillo.h>

//[[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


//--------------------------------------------------------------------------------------------------
//' @title Uniform grid
//' @description 
//' @param a
//' @param b
//' @param N
//' @return Grid.
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
arma::colvec cf_uniform_grid( const double& a, const double& b, const double& N );

//--------------------------------------------------------------------------------------------------
//' @title Exponential grid
//' @description
//' @param l
//' @param a
//' @param b
//' @param n
//' @param E
//' @return Grid.
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
arma::colvec cf_adapt_grid( const double& l, const double& a, const double& b, const double& N,
                            const double& E );

#endif
