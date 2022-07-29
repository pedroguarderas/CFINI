library( testthat )

d <- 100
A <- diag( rep( -2, d ), d )
for ( i in 1:( d - 1 ) ) {
  A[ i, i + 1 ] <- 1
  A[ i + 1, i ] <- 1
}

b <- rep( 1, d )
c <- rep( -1e4, d )
u0 <- rep( 1, d )
n <- 10000
w <- 0.5
e <- 1e-5

S <- cf_psor_solv( u0, A, b, c, w, n, e ) 
u <- S$u

A %*% u - b

us <- solve( A, b )
us

test_that( "Checking numerical solution with PSOR algorithm", {
  expect_lt( norm( us - u, type = '2' ), 0.0001 )
})
