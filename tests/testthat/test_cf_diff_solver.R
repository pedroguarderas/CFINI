library( testthat )

# Time grid
t0 <- 0
t1 <- 1
Nt <- 300
t <- cf_uniform_grid( t0, t1, Nt )

# Space grid
x0 <- -0.5
x1 <- 0.5
Nx <- 150
x <- cf_uniform_grid( x0, x1, Nx )

If <- function( x ) if ( abs(x) <= 0.1 ) return( 1 ) else return( 0 )
If <- Vectorize( If )
I <- sapply( x, FUN = If )
A <- rep( 0, Nt )
B <- rep( 0, Nt )

# Diffusion parameter constant
alpha <- matrix( 10^(-2.3), Nt, Nx )
theta <- 0.25

hx <- ( x1 - x0 ) / ( Nx - 1 )
ht <- ( t1 - t0 ) / ( Nt - 1 )

Ucn <- cf_diff_solv_cns( theta, alpha, I, A, B, t, x )
Ueu <- cf_diff_solv_euls( alpha, I, A, B, t, x )
# persp( t, x, Ucn$u, theta = 80, phi = 45, xlab = 't', ylab = 'x', zlab = 'u', col = 'gold',
#         box = TRUE, axes = TRUE, main = 'Diffusion solution' )

# Computing solution as convolution
alph <- max( alpha )
phi <- function( x ) dnorm( x, mean = 0, sd = sqrt( 2 * alph * t[ Nt ] ) )
s <- sapply( x, FUN = function( y ) integrate( function( x, y ) If( y - x ) * phi( x ), -Inf, Inf, y )$value )
ecn <- matrix( Ucn$u[Nt,] - s, Nx, 1 )
eeu <- matrix( Ueu$u[Nt,] - s, Nx, 1 )

test_that( "Checking Courant-Friedrichs-Lewy condition", {
  expect_lt( ht / ( hx * hx ), 0.5 / max( alpha )  )
})

test_that( "Checking numerical approximation of Crank-Nicolson method", {
  expect_lt( norm( ecn, type = '2' ), 0.5 )
})

test_that( "Checking numerical approximation of Euler implicit method", {
  expect_lt( norm( eeu, type = '2' ), 0.5 )
})

# plot( x, s, type = 'l' )
# points( x, Ueu$u[Nt,], cex = 0.5, pch = 16, col = 'blue' )
# points( x, Ucn$u[Nt,], cex = 0.5, pch = 16, col = 'orange' )
