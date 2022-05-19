library( CFINI )
library( testthat )

n <- 10000
t <- seq( 0, 1, length.out = n )
d <- 2
W <- cf_wiener( d, t )

test_that( "Checking simulation", {
  expect_lte( abs( sqrt( t[2] ) - sd( diff( W[,1] ) ) ), 1e-4 )
})

plot( W[,1], W[,2], type = 'l' )
