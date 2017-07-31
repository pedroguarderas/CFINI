# Pricing options with a multinomial tree ----------------------------------------------------------
#' @title Pricing options with a multinomial tree
#' @description The option pricing can be made with a multinomial tree, this functions can 
#' price different types of options
#' @param Q equivalent discrete martingale measure
#' @param EQ discrete version of the equivalent discrete martingale average
#' @param R term structure of the interest rate, could be a fixed value or a multinomial lattice
#' @param S multinomial lattice
#' @param option function defining the option over S
#' @param type option type a character that specifies the king of option, by default 'E' european
#' option, 'A' american option, 'F' futures option, 'S' swap option, 'P' ...
#' @param option.par list of parameter for the option
#' @return A list with a tree structure of the asset evolution
#' @note Pricing a Forward option is like to price an European option with the identity function
#' regarded as the option
#' @author Pedro Guarderas
#' @importFrom gtools combinations
#' @examples
#' library( gtools )
#' 
#' # Risk neutral average
#' EQ<-function( Q, C ) {
#'   return( sum( Q * C ) )
#' }
#' 
#' # Option
#' call<-function( S, K ) {
#'   max( S - K, 0 )
#' }
#' call.par<-list( K = 99 )
#' 
#' u<-c( 0.9, 1, 1.1 )
#' S0<-100
#' Q<-c( 0.1, 0.7, 0.2 )
#' R<-0.03
#' 
#' S<-MTree( 5, u, S0 )
#' C<-MPricing(  Q, R, EQ, S, call, Type = 'E', call.par )
#' 
#' R<-MTree( 6, u, 0.03 )
#' Cts<-MTreePricing(  Q, R, EQ, S, call, Type = 'E', call.par )
#' @export
CFTreePricing<-function( Q, EQ, R, S, option, Type = 'E', option.par ) {
  C<-S
  N<-length( C )
  n<-length( Q )
  
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
  
  C[[N]]<-sapply( C[[N]], FUN = function( x ) do.call( option, c( x, option.par ) ) )
  
  N<-N-1
  CN<-combinations( n, N, set = TRUE, repeats.allowed = TRUE )
  
  for ( t in N:1 ) {
    CQ<-NULL
    if ( t > 1 )  { 
      
      Ct<-combinations( n, t-1, set = TRUE, repeats.allowed = TRUE )
      
      for ( k in 1:nrow( Ct ) ) {
        J<-Ct[k,]
        K<-NULL
        
        for ( l in 1:n ) {
          r<-sort( c( J, l ) )
          
          for ( i in 1:nrow( CN ) ) {
            if( all( CN[i,] == r ) ) {
              break;
            }
          }
          
          K<-c( K, i )
        }
        
        v<-0
        if ( check == 2 ) {
          v<-1 / ( 1 + R[[t]][k] )
        } else {
          v<-1 / ( 1 + R )
        }
        
        CQ<-c( CQ, EQ( Q, v * C[[t+1]][K] ) )
      }
      
    } else {
      v<-0
      if ( check == 2 ) {
        v<-1 / ( 1 + R[[1]][1] )
      } else {
        v<-1 / ( 1 + R )
      }
      CQ<-EQ( Q, v * C[[2]][1:n] )
    }
    
    if ( Type %in% c( 'E', 'F' ) ) {
      C[[t]]<-CQ
    } else if ( Type == 'A' ) {
      C[[t]]<-sapply( 1:length( CQ ), 
                      FUN = function( k ) max( do.call( option, c( S[[t]][k], option.par ) ), CQ[k] ) )
    }
    
    CN<-Ct
  }
  return( C )
}
