# CFINI --------------------------------------------------------------------------------------------
#' @title CFINI
#' @description This package was created with implement different computational routines that can 
#' be applied to solve pricig problems, that are usually presented in Finance and Actuarial Science.
#' 
#' Here a list of the main numerical methods implemented inside \pkg{CFINI}. Many of this functions 
#' are coded with C++ by employing the packages Rcpp and RcppEigen.
#' \enumerate{
#'   \item Pricing with trees.
#'   \item Pricing with multinomial trees.
#'   \item Ordinary differential equation solver implemented with the predictor-corrector method.
#'   \item Diffusion solver with Euler implicit scheme.
#'   \item Diffusion solver with Crank-Nicolson scheme.
#' }
#' @useDynLib CFINI, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"