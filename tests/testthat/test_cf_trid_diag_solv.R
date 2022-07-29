library( testthat )

d <- 100
A <- diag( rep( -2, d ), d )
for ( i in 1:( d - 1 ) ) {
  A[ i, i + 1 ] <- 1
  A[ i + 1, i ] <- 1
}

b <- rep( 1, d )
c <- rep( -1e3, d )
u0 <- runif( d )
n <- 100
w <- 0.95
e <- 1e-5

S <- cf_psor_solv( u0, A, b, c, w, n, e ) 
u <- S$u

test_that( "Checking numerical solution with PSOR algorithm", {
  expect_lt( norm( A %*% u - b, type = '2' ), 0.0001 )
})