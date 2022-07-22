library( testthat )

# Test predictor-corrector method ------------------------------------------------------------------
n <- 200
t <- seq( 0, 1, length.out = n )
f <- function( t, y ) {
  A <- matrix( c( 1, 0, 0, 1 ), 2, 2 )
  return( A %*% y )
}

v0 <- c( 1, 1 )
f( 1, v0 )
S <- cf_edo_solv_precor( t, v0, f, 5, 1e-6 )

test_that( "Checking numerical solution vs exact solution", {
  expect_lt( norm( exp( rep( 1, 2 ) ) - S$v[ n, 1:2 ], type = '2' ), 0.0001 )
})
