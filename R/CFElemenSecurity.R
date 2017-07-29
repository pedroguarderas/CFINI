# Elementary securities ----------------------------------------------------------------------------
#' @title Elementary securities
#' @description Computes the elementary securities for a multinomial model
#' @param n number of periods for the multinomial model 
#' @param u vector of changes
#' @param r interest rate
#' @param q neutral risk probabilities
#' @return A list with a tree structure of the asset evolution
#' @author Pedro Guarderas
#' @importFrom gtools combinations
#' @examples
#' n<-2
#' q<-c(0.5, 0.5)
#' u<-c(0.9,1.25)
#' r0<-0.03
#' S0<-100
#' P<-ElementarySecurity( n, q, u, r0 )
#' @export
CFElemSecurity<-function( n, q, u, r ) {
  m<-length( q )
  C<-combinations( m, n, set = TRUE, repeats.allowed = TRUE )
  P<-NULL
  for ( i in 1:nrow( C ) ) {
    x<-as.numeric( table( factor( C[i,],  levels = 1:m ) ) )    
    U<-u[ C[i,] ]
    p<-dmultinom( x, prob = rep( 1, m ) ) * ( m^sum(x) ) * prod( q[ C[i,] ] )
    V<-1 / ( 1 + U * r )
    P<-rbind( P, data.table( P = p, U = prod( U ), V = prod( V ) ) )
  }
  
  return( P )
}
