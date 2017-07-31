library( CFINI )

Z<-CFLattice( 4, c(1,1), 100 )
R<-CFLattice( 5, c( 0.9, 1.25 ), 0.06 )

# Zero coupon bond
option<-function( S ) return( S )
EQ<-function( R, Q, C ) {
  return( sum( R * Q * C ) )
}
Q<-c( 0.5, 0.5 )
ZCBlt<-CFLatticePricing( Q, EQ, R, Z, identity, Type = 'E' )


EQ<-function( Q, C ) {
  return( sum( Q * C ) )
}
Z<-CFTree( 4, c(1,1), 100 )
R<-CFTree( 5, c( 0.9, 1.25 ), 0.06 )
ZCBtr<-CFTreePricing( Q, EQ, R, Z, identity, Type = 'E', option.par = list() )

ZCBlt
ZCBtr

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
