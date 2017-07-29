# Multinomial lattice ------------------------------------------------------------------------------
#' @title Multinomial lattice
#' @description The option pricing can be made with a multinomial lattice, this functions can 
#' generate such a lattice with aid of some give parameters.
#' @param N number of periods for the multinomial model 
#' @param U vector of changes
#' @param S0 initial price
#' @return A list with a tree structure of the asset evolution
#' @author Pedro Guarderas
#' @examples
#' s<-0.3
#' T<-0.25
#' N<-15
#' u<-exp( s * sqrt( T / N ) )
#' d<-1/u
#' r<-0.02
#' S0<-100
#' R<-exp( r * T / N )
#' q<-( R - d ) / ( u - d )
#' U<-c( d, u )
#' S<-MLattice( N, U, S0 )
#' @export
CFLattice<-function( N, U, S0 ) {
  S<-list( S0, U * S0 )
  n<-length( U )
  
  for ( t in 3:(N+1) ) {
    M<-n + (t-2)*(n-1)
    K<-M - n + 1
    S[[t]]<-c( U * S[[t-1]][1], U[n] * S[[t-1]][2:K] )
  }
  return( S )
}
