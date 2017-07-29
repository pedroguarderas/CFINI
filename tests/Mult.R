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


