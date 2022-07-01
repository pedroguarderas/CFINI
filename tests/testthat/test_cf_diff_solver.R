library( testthat )

# Time grid
t0 <- 0
t1 <- 1
Nt <- 200
t <- cf_uniform_grid( t0, t1, Nt )

# Space grid
x0 <- -0.5
x1 <- 0.5
Nx <- 200
x <- cf_uniform_grid( x0, x1, Nx )

If <- function( x ) if ( abs(x) <= 0.1 ) return( 1 ) else return( 0 )
If <- Vectorize( If )
I <- sapply( x, FUN = If )
A <- rep( 0, Nt )
B <- rep( 0, Nt )

alpha <- 10^(-2.3)
theta <- 0.5

hx <- ( x1 - x0 ) / ( Nx - 1 )
ht <- ( t1 - t0 ) / ( Nt - 1 )

# print( ht / ( hx * hx ) )
# print( 0.5 / alpha )

U <- cf_diff_solv_cns( alpha, theta, I, A, B, t, x )
# persp( t, x, U$u, theta = 80, phi = 45, xlab = 't', ylab = 'x', zlab = 'u', col = 'gold',
#        box = TRUE, axes = TRUE, main = 'Diffusion solution' )

phi <- function( x ) dnorm( x, mean = 0, sd = sqrt( 2 * alpha * t[ Nt ] ) )

# Computing solution as convolution
s <- sapply( x, FUN = function( y ) integrate( function( x, y ) If( y - x ) * phi( x ), -Inf, Inf, y )$value )
e <- matrix( U$u[Nt,] - s, Nx, 1 )

test_that( "Verification of numerical solution with Crank-Nicolson method", {
  expect_lt( norm( e, type = '2' ), 0.5 )
})
