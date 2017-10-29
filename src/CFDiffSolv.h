#ifndef __CF_DIFF_SOLV__
#define __CF_DIFF_SOLV__

#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include "CFTriDiagSolv.h"

using namespace Rcpp;

//--------------------------------------------------------------------------------------------------
//' @title Diffusion solver
//' @description Solver for standard diffusion problems
//' @param alpha
//' @param I
//' @param A
//' @param B
//' @param t
//' @param x
//' @return List with solution parameters
//' @note Diffusion solver for pricing options
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
List CFDiffSolvES( const double& alpha,
                   const arma::colvec& I,
                   const arma::colvec& A,
                   const arma::colvec& B,
                   const arma::colvec& t,
                   const arma::colvec& x );

//--------------------------------------------------------------------------------------------------
//' @title Diffusion solver
//' @description Solver for diffusion problems implemented with Crank-Nicolson scheme
//' @param alpha
//' @param theta
//' @param I
//' @param A
//' @param B
//' @param t
//' @param x
//' @return List with solution parameters
//' @note Diffusion solver for pricing options
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
List CFDiffSolvCNS( const double& alpha,
                    const double& theta,
                    const arma::colvec& I,
                    const arma::colvec& A,
                    const arma::colvec& B,
                    const arma::colvec& t,
                    const arma::colvec& x );

//--------------------------------------------------------------------------------------------------
//' @title Black-Scholes solver
//' @description Solver for Black-Scholes models implemented with Crank-Nicolson scheme
//' @param sigma
//' @param rate
//' @param theta
//' @param I
//' @param A
//' @param B
//' @param t
//' @param x
//' @return List with solution parameters
//' @note pricing options
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
List CFBlackScholesSolvCNS( const double& sigma,
                            const double& rate,
                            const double& theta,
                            const arma::colvec& I,
                            const arma::colvec& A,
                            const arma::colvec& B,
                            const arma::colvec& t,
                            const arma::colvec& x );

#endif
