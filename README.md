
<!-- README.md is generated from README.Rmd. Please edit that file -->

\#CFINI

This package was created with the pourpose to implement different
computational routines that can be applied to solve stochastic problems,
mainly with applications in finance and actuarial science.

Several functions are implemented with C++ by employing the packages
Rcpp and RcppEigen.

-   Pricing with trees
-   Pricing with multinomial trees
-   Ordinary differential equation solver implemented with the
    predictor-corrector method
-   Diffusion solver with Euler implicit scheme
-   Diffusion solver with Crank-Nicolson scheme

The different methods are implemented following the following sources of
information:

-   The book [Finite difference methods in financial
    engineering](https://www.wiley.com/en-us/Finite+Difference+Methods+in+Financial+Engineering:+A+Partial+Differential+Equati%20on+Approach-p-9781118856482)
    of Daniel J. Duffy.
-   The specialization of coursera [Financial Engineering and Risk
    Management](https://www.coursera.org/specializations/financialengineering)
    a MOOC provided by Columbia University.
-   The MOOC of coursera [Interest Rate
    Models](https://www.coursera.org/learn/interest-rate-models)
    provided by EPFL.
