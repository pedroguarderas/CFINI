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
//' @param alpha matrix discretization Diffusion parameter
//' @param u0 vector with discretization of the initial condition
//' @param u1 vector with discretization of the inferior boundary condition
//' @param u2 vector with discretization of the superior boundary condition
//' @param t vector with discretization of the time grid
//' @param x vector with discretization of the spatial grid
//' @param is_initial specifies if the condition is final of initial
//' @return List with solution parameters
//' @note Diffusion solver for pricing options
//' @author Pedro Guarderas
//' \email{pedro.felipe.guarderas@@gmail.com}
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
//' @param alpha matrix discretization Diffusion parameter
//' @param u0 vector with discretization of the initial condition
//' @param u1 vector with discretization of the inferior boundary condition
//' @param u2 vector with discretization of the superior boundary condition
//' @param t vector with discretization of the time grid
//' @param x vector with discretization of the spatial grid
//' @param is_initial specifies if the condition is final of initial
//' @return List with solution parameters
//' @note Diffusion solver for pricing options
//' @author Pedro Guarderas
//' \email{pedro.felipe.guarderas@@gmail.com}
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
//' @param u0 vector with discretization of the initial condition
//' @param u1 vector with discretization of the inferior boundary condition
//' @param u2 vector with discretization of the superior boundary condition
//' @param t vector with discretization of the time grid
//' @param x vector with discretization of the spatial grid
//' @return List with solution parameters
//' @note pricing options
//' @author Pedro Guarderas
//' \email{pedro.felipe.guarderas@@gmail.com}
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
