# Multinomial tree ---------------------------------------------------------------------------------
#' @title Multinomial tree
#' @description The option pricing can be made with a multinomial tree, this functions can 
#' generate such a lattice with aid of some give parameters.
#' @param n number of periods for the multinomial model 
#' @param u vector of changes
#' @param S0 initial price
#' @return A list with a tree structure of the asset evolution
#' @author Pedro Guarderas
#' \email{pedro.felipe.guarderas@@gmail.com}
#' @importFrom gtools combinations
#' @export
cf_tree <- function( n, u, S0 ) {
  m <- length( u )
  
  tree <- list( S0, S0 * u )
  if ( n >= 2 ) {
    for ( k in 2:n ) {
      # Generator of combinations with replacement
      c <- rep( 1, k )
      j <- k
      p <- prod( u[ c ] ) * S0
      while ( j > 0 ) {
        l <- k
        while ( l >= j ) {
          while ( c[ l ] < m ) {
            c[ l ] <- c[ l ] + 1
            p <- c( p, prod( u[ c ] ) * S0 )
          }
          l <- l - 1
          if ( l > 0 ) {
            if ( c[l] < m ) {
              c[ l:k ] <- c[ l ] + 1
              p <- c( p, prod( u[ c ] ) * S0 )
              l <- k
            }
          }
        } 
        j <- j - 1 
      }
      
      tree[[ k + 1 ]] <- p
    }
  }
  return( new( 'cflattice', lattice = tree ) )
}
