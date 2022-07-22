library( CFINI )
library( testthat )

n <- 1e4
t <- seq( 0, 1, length.out = n )
d <- 2
W <- cf_wiener( d, t )

test_that( "Checking simulation", {
  expect_lte( abs( sqrt( t[2] ) - sd( diff( W[,1] ) ) ), 1e-4 )
})

plot( W[,1], W[,2], type = 'l' )

# Test predictor-corrector method ------------------------------------------------------------------
n <- 1000
t <- seq( 0, 1, length.out = n )
f <- function( t, y ) {
 A <- matrix( c( 1, 0, 0, 1 ), 2, 2 )
 return( A %*% y )
}

v0 <- c( 1, 1 )
f( 1, v0 )
S <- cf_edo_solv_precor( t, v0, f, 10, 1e-6 )
exp( 1 )
S$v[ n ]
