#ifndef __CF_DIFF_SOLV__
#define __CF_DIFF_SOLV__

#include <RcppEigen.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

#include "cf_tri_diag_solv.h"

using namespace Rcpp;

//--------------------------------------------------------------------------------------------------
//' @title Diffusion solver with implicit Euler scheme
//' @description Solver for diffusion problems implemented using the Euler implicit scheme
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
List cf_diff_solv_euls( const Eigen::MatrixXd& alpha,
                        const Eigen::VectorXd& u0,
                        const Eigen::VectorXd& u1,
                        const Eigen::VectorXd& u2,
                        const Eigen::VectorXd& t,
                        const Eigen::VectorXd& x,
                        const bool is_initial );

//--------------------------------------------------------------------------------------------------
//' @title Diffusion solver with Crank-Nicolson scheme
//' @description Solver for diffusion problems implemented with Crank-Nicolson scheme
//' @param theta Parameter for the Crank-Nicolson scheme
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
List cf_diff_solv_cns( const double& theta,
                       const Eigen::MatrixXd& alpha,
                       const Eigen::VectorXd& u0,
                       const Eigen::VectorXd& u1,
                       const Eigen::VectorXd& u2,
                       const Eigen::VectorXd& t,
                       const Eigen::VectorXd& x,
                       const bool is_initial );

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
List cf_black_scholes_solv_cns( const double& sigma,
                                const double& rate,
                                const double& theta,
                                const Eigen::VectorXd& u0,
                                const Eigen::VectorXd& u1,
                                const Eigen::VectorXd& u2,
                                const Eigen::VectorXd& t,
                                const Eigen::VectorXd& x );

#endif
