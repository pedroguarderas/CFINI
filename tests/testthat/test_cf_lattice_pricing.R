library( testthat )

s <- 0.3
T <- 0.25
N <- 15
c <- 0.01
u <- exp( s * sqrt( T / N ) )
d <- 1 / u
r <- 0.02
S0 <- 100
K <- 110
R <- exp( r * T / N )
q <- ( R - d ) / ( u - d )
Q <- c( 1 - q, q )
U <- c( d, u )

# K is global
call <- function( S ) {
  max( S - K, 0 )
}

put <- function( S ) {
  max( K - S, 0 )
}

# Equivalent measure
EQ <- function( R, Q, C ) {
  return( sum( R * Q * C ) )
}

S <- cf_lattice( N, U, S0 )

# Pricing zero coupon bond
ZCB <- cf_lattice_pricing( Q, EQ, R, S, identity, type = 'A' )
test_that( "Checking expansion for lattice pricing", {
  expect_equal( abs( ZCB[[16]][16] - U[2]^15 * S0 ), 0.0 )
})

# Pricing american call
Ca <- cf_lattice_pricing( Q, EQ, R, S, call, type = 'A' )

# Pricing american put
Pa <- cf_lattice_pricing( Q, EQ, R, S, put, type = 'A' )
