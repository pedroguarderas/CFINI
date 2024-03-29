# Computing amortization ---------------------------------------------------------------------------
#' @title Computing amortization
#' @description Function to compute the amortization tables
#' @param r interest rate
#' @param y number of years
#' @param p number of periods for year
#' @param M loan ammount
#' @return A data.table with amortization information
#' @author Pedro Guarderas
#' @examples
#' r <- 0.16
#' y <- 10
#' p <- 12
#' M <- 10000
#' A <- cf_amortization( r, y, p, M )
#' sum( A$c )
#' sum( A$i )
#' @importFrom data.table := data.table
#' @export
cf_amortization <- function( r, y, p, M ) {
  
  I <- r / p
  N <- y * p
  R <- -( 1 - ( 1 + I )^( N + 1 ) ) / I
  C <- M * ( I * ( 1 + I )^N ) / ( ( 1 + I )^N - 1)
  
  t <- seq( 0, N - 1, 1 )
  m <- M * ( 1 + I )^t  + C * ( 1 - ( 1 + I )^t ) / I
  A <- data.table( 
    t = t, 
    y = t / p,
    r = ( 1 + I )^t,
    v = ( 1 + I )^(-t),
    m = m,
    i = m * I,
    c = C - m * I )
  
  return( A )
  
}
