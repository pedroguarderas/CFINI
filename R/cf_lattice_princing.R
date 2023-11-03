# Option pricing with the multinomial model --------------------------------------------------------
#' @title Option pricing with the multinomial model
#' @description Function for pricing an option with a multinomial model
#' @param Q equivalent discrete martingale measure
#' @param EQ discrete version of the equivalent discrete martingale average
#' @param R term structure of the interest rate, could be a fixed value or a multinomial lattice
#' @param S multinomial lattice
#' @param option function defining the option over S
#' @param type option type a character that specifies the king of option, by default 'E' european
#' option, 'A' american option, 'F' futures option, 'S' swap option, 'P' ...
#' @return A list with a tree structure of the asset evolution
#' @author Pedro Guarderas
#' @seealso \code{\link{cflattice-class}}
#' @examples
#' s <- 0.3
#' T <- 0.25
#' N <- 15
#' c <- 0.01
#' u <- exp( s * sqrt( T / N ) )
#' d <- 1 / u
#' r <- 0.02
#' S0 <- 100
#' K <- 110
#' R <- exp( r * T / N )
#' q <- ( R - d ) / ( u - d )
#' Q <- c( 1 - q, q )
#' U <- c( d, u )
#' 
#' # K is global
#' call <- function( S ) {
#'   max( S - K, 0 )
#' }
#' 
#' put <- function( S ) {
#'   max( K - S, 0 )
#' }
#' 
#' # Equivalent measure
#' EQ <- function( R, Q, C ) {
#'   return( sum( R * Q * C ) )
#' }
#' 
#' S <- cf_lattice( N, U, S0 )
#' 
#' # Pricing zero coupon bond
#' ZCB <- cf_lattice_pricing( Q, EQ, R, S, identity, type = 'A' )
#' 
#' # Pricing american call
#' Ca <- cf_lattice_pricing( Q, EQ, R, S, call, type = 'A' )
#' 
#' # Pricing american put
#' Pa <- cf_lattice_pricing( Q, EQ, R, S, put, type = 'A' )
#' @export
cf_lattice_pricing <- function( Q, EQ, R, S, option, type = 'E' ) {
  C <- S
  N <- length(S)
  n <- length(Q)
  
  check <- 1
  if ( length( R ) > 1 ) {
    ns <- unlist( lapply( S, FUN = length ) )
    nr <- unlist( lapply( R, FUN = length ) )
    if ( length( ns ) < length(nr) ) {
      if ( all( ns == nr[1:length(ns)] ) ) check <- 2
    } else if ( length( ns ) == length(nr) ) {
      if ( all( ns == nr ) ) check <- 2
    } else if ( length( ns ) > ( length(nr) + 1 ) ) {
      if ( all( ns[1:length(nr)] == nr ) ) check <- 2
    }
  }
  
  if ( type == 'E' & check == 1 ) {
    C@lattice[[ N ]] <- sapply( C@lattice[[ N ]], FUN = option )
    for ( t in (N-1):1 ) {
      M <- n + (t-2)*(n-1)
      C@l[[ t ]] <- sapply( 
        1:M, 
        FUN = function(k) EQ( 1/R, Q, C@lattice[[ t + 1 ]][k:(k+n-1)] )  )
    }
  } else if ( type == 'E' & check == 2 ) {
    C@lattice[[ N ]] <- sapply( C@lattice[[ N ]], FUN = option )
    for ( t in (N-1):1 ) {
      M <- n + (t-2)*(n-1)
      C@lattice[[ t ]] <- sapply( 
        1:M, 
        FUN = function(k) EQ( 1 / ( 1 + R@lattice[[ t ]][k] ), Q, C@lattice[[ t + 1 ]][k:(k+n-1)] )  )
    }
  } else if ( type == 'A' & check == 1 ) {
    C@lattice[[ N ]] <- sapply( C@lattice[[ N ]], FUN = option )
    for ( t in (N-1):1 ) {
      M <- n + (t-2)*(n-1)
      C@lattice[[ t ]] <- sapply( 
        1:M, 
        FUN = function(k) max( option( S@lattice[[ t ]][k] ),EQ( 1/R, Q, C@lattice[[ t + 1 ]][k:(k+n-1)] ) ) )
    }
  } else if ( type == 'A' & check == 2 ) {
    C@lattice[[ N ]] <- sapply( C@lattice[[ N ]], FUN = option )
    for ( t in (N-1):1 ) {
      M <- n + (t-2)*(n-1)
      C@lattice[[ t ]] <- sapply( 
        1:M, 
        FUN = function(k) max( option( S@lattice[[ t ]][k] ), EQ( 1 / ( 1 + R@lattice[[ t ]][k] ), Q, C@lattice[[ t + 1 ]][k:(k+n-1)] ) ) )
    }
  } else if ( type == 'F' & check == 1 ) {
    for ( t in (N-1):1 ) {
      M <- n + (t-2)*(n-1)
      C@lattice[[ t ]] <- sapply( 
        1:M, 
        FUN = function(k) EQ( 1 / R, Q, C@lattice[[ t + 1 ]][ k:(k+n-1) ] )  )
    }
  } else if ( type == 'F' & check == 2 ) {
    for ( t in (N-1):1 ) {
      M <- n + (t-2)*(n-1)
      C@lattice[[ t ]] <- sapply( 
        1:M, 
        FUN = function(k) EQ( 1 / ( 1 + R@lattice[[ t ]][k] ), Q, C@lattice[[ t + 1 ]][ k:(k+n-1) ] )  )
    }
  } else if ( type == 'S' & check == 2 ) {
    C@lattice[[ N ]] <- sapply( C@lattice[[ N ]], FUN = option )
    C@lattice[[ N ]] <- C@lattice[[ N ]] / ( 1 + R[N] )
    for ( t in (N-1):1 ) {
      M <- n + (t-2)*(n-1)
      C@lattice[[ t ]] <- sapply( 1:M, FUN = function(k) {
        EQ( 1 / ( 1 + R@lattice[[ t ]][k] ), Q, C@lattice[[ t + 1 ]][k:(k+n-1)] ) +
          option( S@lattice[[ t ]][k] ) / ( 1 + R@lattice[[ t ]][k] )
      } )
    }
  } # else if ( type = 'P' & check == 2 ) {
  #   C[1] <- 1
  #   for ( t in 2:N ) {
  #     M <- n + (t-2)*(n-1)
  #     C@lattice[[ t ]] <- sapply( 1:M, FUN = function(k) {
  #       EQ( 1 / ( 1 + R@lattice[[ t ]][k] ), Q[], C[t-1][k:(k+n-1)] )
  #       
  #     } )
  #   }
  #   
  # }
  return( C )
}
