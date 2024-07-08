
<!-- badges: start -->

[![R-CMD-check](https://github.com/pedroguarderas/CFINI/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pedroguarderas/CFINI/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

# CFINI

## Introduction

The CFINI package implements routines for the solution of pricing problems usually encountered in
Financial and Actuarial Science.

<!-- ?? include? are there more? -->
The numerical methods implemented in CFINI include: 

<!-- ?? Check names of methods below -->
* Pricing with trees

* Pricing with multinomial trees

* Ordinary differential equation solver implemented with the predictor-corrector method

* Diffusion solver with implicit Euler method

* Diffusion solver with Crank-Nicolson method

Many functions are coded with C++ by employing the [Rcpp][] and [RcppEigen][] packages.

[Rcpp]: https://cran.r-project.org/web/packages/Rcpp/index.html

[RcppEigen]: https://cran.r-project.org/web/packages/RcppEigen/index.html

### References

* Shreve, Steven E., *Stochastic Calculus for Finance I: The Binomial Asset Pricing Model*
  (1st. edn, [Springer][springer_1]).

* Duffy, Daniel J., *Finite Difference Methods in Financial Engineering: A Partial Differential
  Equation Approach* (2013, [Wiley][wiley_1]).

* Haugh, Hirsa, Iyengar, *Columbia University MOOC: Financial Engineering and Risk Management
  Specialization* (2024,[Coursera][coursera_1]).

* Filipović D., *École Polytechnique Fédérale de Lausanne MOOC: Interest Rate Models* (2024,
  [Coursera][coursera_2]).

[springer_1]: https://doi.org/10.1007/978-0-387-22527-2

[wiley_1]: https://www.wiley.com/en-us/Finite+Difference+Methods+in+Financial+Engineering:+A+Partial+Differential+Equati%20on+Approach-p-9781118856482

[coursera_1]: https://www.coursera.org/specializations/financialengineering

[coursera_2]: https://www.coursera.org/learn/interest-rate-models

## Installation

To install **CFINI** directly from GitHub use `devtools`:

```R

library( devtools )
install_github( "pedroguarderas/CFINI" )

```
