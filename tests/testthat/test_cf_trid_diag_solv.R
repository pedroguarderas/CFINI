library( testthat )

d <- 50
A <- diag( rep( -2, d ), d )
for ( i in 1:( d - 1 ) ) {
  A[ i, i + 1 ] <- 1
  A[ i + 1, i ] <- 1
}

b <- rep( 1, d )
c <- rep( -1e4, d )
u0 <- runif( d )
n <- 100000
w <- 0.8
e <- 1e-6

S <- cf_psor_solv( u0, A, b, c, w, n, e ) 
u <- S$u

test_that( "Checking numerical solution with PSOR algorithm", {
  expect_lt( norm( A %*% u - b, type = '2' ), 1e-6 )
})
