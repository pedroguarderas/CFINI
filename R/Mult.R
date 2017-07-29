library( CFINI )
Z<-MLattice( 4, c(1,1), 100 )
R<-MLattice( 2, c( 0.9, 0.99, 1.1 ), 0.05 )

# Zero coupon bond
identity<-function( S ) return( S )
EQ<-function(R,Q,C) {
  return( sum( R * Q * C ) )
}
Q<-c( 0.5, 0.5 )
ZCB<-MPricing( Z, identity, EQ, R, Q, Type = 'E' )

# Put call
K<-84
call<-function( S ) max( S - K, 0 )
n<-3
n<-n+1
CZCB<-MPricing( ZCB[1:n], call, EQ, R, Q, Type = 'E' )

# Put call
K<-88
put<-function( S ) max( K - S, 0 )
n<-3
n<-n+1
PZCB<-MPricing( ZCB[1:n], put, EQ, R, Q, Type = 'A' )

# Caplet
K<-0.02
n<-4
n<-n+1
R<-MLattice( 5, c( 0.9, 1.25 ), 0.06 )
Z<-lapply( R, FUN = function( r ) sapply( r, FUN = function( S ) call( S ) / ( 1 + S ) ) )
CR<-MPricing( Z, identity, EQ, R, Q, Type = 'E' )


n<-2
combn( 1:3, 2 )
expand.grid( 1:n, 1:n )
q<-c(0.2,0.5,0.3)
x<-c(1,1,2)
dmultinom( x, prob = rep( 1, length(x))/length(x) ) * (length(x)^sum(x))
q^x


library( gtools )
library( data.table )

MTree<-function( n, u, S0 ) {
  m<-length( u )
  N<-n
  
  Tree<-list( S0 )
  for ( k in 1:N ) {
    C<-combinations( m, k, set = TRUE, repeats.allowed = TRUE )
    U<-NULL
    for ( i in 1:nrow( C ) ) {
      x<-as.numeric( table( factor( C[i,],  levels = 1:m ) ) )    
      U<-c( U, prod( u[ C[i,] ] ) )
    }  
    Tree[[k+1]]<-S0 * U
  }
  return( Tree )
}

ElementarySecurity<-function( n, q, u, r ) {
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

n<-2
q<-c(0.5, 0.5)
u<-c(0.9,1.25)
r0<-0.03
S0<-100
P<-ElementarySecurity( n, q, u, r0 )
P[ , S := S0 * P * V ]
sum( P$S )
Tree<-MTree( 1, c(1,1), S0 )

R<-MTree( 4, u, 0.06 )
MPricing( MTree( 4, c(1,1), S0 ), identity, EQ, R, q, Type = 'E' )
