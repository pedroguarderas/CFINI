# Option pricing with the multinomial model --------------------------------------------------------
#' @title Option pricing with the multinomial model
#' @description Function for pricing an option with a multinomial model
#' @param S multinomial lattice
#' @param option function defining the option over S
#' @param EQ discrete version of the equivalent discrete martingale average
#' @param R term structure of the interest rate, could be a fixed value or a multinomial lattice
#' @param Q equivalent discrete martingale measure
#' @param type option type a character that specifies the king of option, by default 'E' european
#' option, 'A' american option, 'F' futures option, 'S' swap option, 'P' ...
#' @return A list with a tree structure of the asset evolution
#' @author Pedro Guarderas
#' @seealso \code{\link{MLattice}}
#' @examples
#' s<-0.3
#' T<-0.25
#' N<-15
#' c<-0.01
#' u<-exp( s * sqrt( T / N ) )
#' d<-1/u
#' r<-0.02
#' S0<-100
#' K<-110
#' R<-exp( r * T / N )
#' q<-( R - d ) / ( u - d )
#' Q<-c( 1-q, q )
#' U<-c( d, u )
#' 
#' # K is global
#' call<-function( S ) {
#'   max( S - K, 0 )
#' }
#' 
#' put<-function( S ) {
#'   max( K - S, 0 )
#' }
#'
#' # Equivalent measure
#' EQ<-function(R,Q,C) {
#'   return( sum( R * Q * C ) )
#' }
#' 
#' S<-CFLattice( N, U, S0 )
#' 
#' # Pricing american call
#' Ca<-CFLatticePricing( S, call, EQ, R, Q, Type = 'A' )
#' 
#' # Pricing american put
#' Pa<-CFLatticePricing( S, put, EQ, R, Q, Type = 'A' )
#' @export
CFLatticePricing<-function( S, option, EQ, R, Q, Type = 'E' ) {
  C<-S
  N<-length(S)
  n<-length(Q)
  
  check<-1
  if ( length( R ) > 1 ) {
    ns<-unlist( lapply( S, FUN = length ) )
    nr<-unlist( lapply( R, FUN = length ) )
    if ( length( ns ) < length(nr) ) {
      if ( all( ns == nr[1:length(ns)] ) ) check<-2
    } else if ( length( ns ) == length(nr) ) {
      if ( all( ns == nr ) ) check<-2
    } else if ( length( ns ) > ( length(nr) + 1 ) ) {
      if ( all( ns[1:length(nr)] == nr ) ) check<-2
    }
  }
  
  if ( Type == 'E' & check == 1 ) {
    C[[N]]<-sapply( C[[N]], FUN = option )
    for ( t in (N-1):1 ) {
      M<-n + (t-2)*(n-1)
      C[[t]]<-sapply( 1:M, FUN = function(k) EQ( 1/R, Q, C[[t+1]][k:(k+n-1)] )  )
    }
  } else if ( Type == 'E' & check == 2 ) {
    C[[N]]<-sapply( C[[N]], FUN = option )
    for ( t in (N-1):1 ) {
      M<-n + (t-2)*(n-1)
      C[[t]]<-sapply( 1:M, FUN = function(k) EQ( 1 / ( 1 + R[[t]][k] ), Q, C[[t+1]][k:(k+n-1)] )  )
    }
  } else if ( Type == 'A' & check == 1 ) {
    C[[N]]<-sapply( C[[N]], FUN = option )
    for ( t in (N-1):1 ) {
      M<-n + (t-2)*(n-1)
      C[[t]]<-sapply( 1:M, FUN = function(k) max( option( S[[t]][k] ),
                                                  EQ( 1/R, Q, C[[t+1]][k:(k+n-1)] ) ) )
    }
  } else if ( Type == 'A' & check == 2 ) {
    C[[N]]<-sapply( C[[N]], FUN = option )
    for ( t in (N-1):1 ) {
      M<-n + (t-2)*(n-1)
      C[[t]]<-sapply( 1:M, FUN = function(k) max( option( S[[t]][k] ),
                                                  EQ( 1 / ( 1 + R[[t]][k] ), Q, C[[t+1]][k:(k+n-1)] ) ) )
    }
  } else if ( Type == 'F' & check == 1 ) {
    for ( t in (N-1):1 ) {
      M<-n + (t-2)*(n-1)
      C[[t]]<-sapply( 1:M, FUN = function(k) EQ( 1 / R, Q, C[[t+1]][k:(k+n-1)] )  )
    }
  } else if ( Type == 'F' & check == 2 ) {
    for ( t in (N-1):1 ) {
      M<-n + (t-2)*(n-1)
      C[[t]]<-sapply( 1:M, FUN = function(k) EQ( 1 / ( 1 + R[[t]][k] ), Q, C[[t+1]][k:(k+n-1)] )  )
    }
  } else if ( Type == 'S' & check == 2 ) {
    C[[N]]<-sapply( C[[N]], FUN = option )
    C[[N]]<-C[[N]] / ( 1 + R[[N]] )
    for ( t in (N-1):1 ) {
      M<-n + (t-2)*(n-1)
      C[[t]]<-sapply( 1:M, FUN = function(k) {
        EQ( 1 / ( 1 + R[[t]][k] ), Q, C[[t+1]][k:(k+n-1)] ) +
          option( S[[t]][k] ) / ( 1 + R[[t]][k] )
      } )
    }
  } # else if ( Type = 'P' & check == 2 ) {
  #   C[[1]]<-1
  #   for ( t in 2:N ) {
  #     M<-n + (t-2)*(n-1)
  #     C[[t]]<-sapply( 1:M, FUN = function(k) {
  #       EQ( 1 / ( 1 + R[[t]][k] ), Q[], C[[t-1]][k:(k+n-1)] )
  #       
  #     } )
  #   }
  #   
  # }
  return( C )
}
