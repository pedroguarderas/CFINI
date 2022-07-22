#ifndef __CF_ODE_SOLVER__
#define __CF_ODE_SOLVER__

#include <RcppEigen.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

//--------------------------------------------------------------------------------------------------
//' @title Ordinary differenatial equation solver
//' @description Solver for ordinary differential equations, implemented with the predictor 
//' corrector method.
//' @param t time grid, could be adapted
//' @param v0 initial value
//' @param f dynamic function
//' @return Return the numerical solution to the ODE.
//' @note The solver is implemented to compute the solution forwards in time.
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
List cf_edo_solv_precor( const Eigen::VectorXd& t,
                         const Eigen::VectorXd& v0,
                         Function f,
                         const int& m = 2,
                         const double& err = 1e-2 );

#endif