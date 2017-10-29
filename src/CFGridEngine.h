#ifndef __CF_GRID_ENGINE__
#define __CF_GRID_ENGINE__

#include <RcppArmadillo.h>

//[[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


//--------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
arma::colvec GridUniform( const double& a, const double& b, const double& N );

//--------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
arma::colvec GridExpAddapt( const double& l, const double& a, const double& b, const double& N,
                            const double& E );

#endif
