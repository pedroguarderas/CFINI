library( testthat )

d <- 100
A <- diag( rep( -2, d ), d )
for ( i in 1:( d - 1 ) ) {
  A[ i, i + 1 ] <- 1
  A[ i + 1, i ] <- 1
}

b <- rep( 1, d )
c <- rep( -1e6, d )
u0 <- rep( 1, d )
n <- 1000
w <- 0.5
e <- 1e-4

S <- cf_sor_solv( u0, A, b, c, w, n, e ) 
u <- S$u

A %*% u - b
