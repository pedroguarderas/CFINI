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
//' @param alpha Diffusion parameter
//' @param I Initial condition
//' @param A Inferior boundary condition
//' @param B Superior boundary condition
//' @param t Time grid
//' @param x Spatial grid
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
//' @title Diffusion solver with Crank-Nicolson scheme
//' @description Solver for diffusion problems implemented with Crank-Nicolson scheme
//' @param alpha Diffusion parameter
//' @param theta Parameter for the Crank-Nicolson scheme
//' @param I Initial condition
//' @param A Inferior boundary condition
//' @param B Superior boundary condition
//' @param t Time grid
//' @param x Spatial grid
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
//' @description Solver for Black-Scholes models implemented with the Crank-Nicolson numerical 
//' scheme
//' @param sigma Volatility
//' @param rate Interest rate
//' @param theta Parameter for the Crank-Nicolson scheme
//' @param I Initial condition
//' @param A Inferior boundary condition
//' @param B Superior boundary condition
//' @param t Time grid
//' @param x Spatial grid
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
