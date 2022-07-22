# Multinomial tree ---------------------------------------------------------------------------------
#' @title Multinomial tree
#' @description The option pricing can be made with a multinomial tree, this functions can 
#' generate such a lattice with aid of some give parameters.
#' @param n number of periods for the multinomial model 
#' @param u vector of changes
#' @param S0 initial price
#' @return A list with a tree structure of the asset evolution
#' @author Pedro Guarderas
#' @importFrom gtools combinations
#' @export
cf_tree <- function( n, u, S0 ) {
  m <- length( u )
  N <- n
  
  Tree <- list( S0 )
  for ( k in 1:N ) {
    C <- combinations( m, k, set = TRUE, repeats.allowed = TRUE )
    U <- NULL
    for ( i in 1:nrow( C ) ) {
      U <- c( U, prod( u[ C[i,] ] ) )
    }  
    Tree[[k+1]] <- S0 * U
  }
  return( Tree )
}