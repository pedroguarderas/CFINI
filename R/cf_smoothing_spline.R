# Spline basado en el método de Lorimier -----------------------------------------------------------
#' @title Spline basado en el método de Lorimier
#' @description Spline basado en el método de Lorimier
#' @param t tiempo
#' @param y rendimiento
#' @param alpha factor de ajuste
#' @return lista con resultados de la interpolación, coeficientes y matriz
#' @author Pedro Guarderas
#' @import Matrix
#' @export
cf_smoothing_spline<-function( t, y, alpha ) {

  N<-length( t ) + 1
  A<-Matrix( 0, N, N )
  A[1,]<-c( 0, t )
  A[,1]<-c( 0, t )

  IP<-function( s, t ) {
    m<-min( s, t )
    return( s * t + s * t * m - ( s + t ) * ( m^2 ) / 2 + ( m^3 ) / 3 )
  }

  for ( i in 1:(N-1) ) {
    A[i+1,2:N]<-sapply( t, FUN = IP, t = t[i] )
  }

  A[2:N,1:N]<-alpha * A[2:N,1:N]
  diag( A )<-diag( A ) + 1
  A[1,1]<-0

  u<-c( 0, alpha * y * t )
  b<-solve( A, u )

  return( list( b = b, A = A ) )
}

# Tasa ---------------------------------------------------------------------------------------------
#' @title Tasa
#' @description Tasa
#' @param u valor
#' @param b coeficientes de interpolación
#' @param t tiempo
#' @return vector
#' @author Pedro Guarderas
#' @export
cf_forward_rate<-function( u, t, b ) {
  spline_base<-function( s, t ) {
    m<-min( s, t )
    return( t * m - 0.5 * m^2 + t )
  }

  N<-length( t )
  f<-matrix( 0, length( u ), length( b ) )
  f[,1]<-1
  for ( i in 1:N ) {
    f[,i+1]<-sapply( u, FUN = spline_base, t = t[i] )
  }

  f<-as.numeric( f %*% b )
  return(f)
}

# Rendimiento --------------------------------------------------------------------------------------
#' @title Rendimiento
#' @description Rendimiento
#' @param u valor
#' @param b coeficientes de interpolación
#' @param t tiempo
#' @return vector
#' @author Pedro Guarderas
#' @export
cf_yield<-function( u, b, t ) {
  # Integration base function
  spline_base_integration<-function( t, u ) {
    A<-0
    B<-0
    if ( u <= t ) {
      A<-0.5 * u^2
      B<-(1/3) * u^3
    } else {
      A<-0.5 * t^2 + ( u - t ) * t
      B<-(1/3) * t^3 + ( u - t ) * t^2
    }
    I<- t * A - 0.5 * B + t * u
    return( I )
  }

  n<-length( b )
  y<-sapply( t, FUN = spline_base_integration, u = u )
  y<-b[1] + sum( b[2:n] * y ) / u
  return( y )
}
