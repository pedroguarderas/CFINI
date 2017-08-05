library( Rcpp )
library( CFINI )

PointMod<-Module( 'PointMod', PACKAGE = 'CFINI' )
Point<-PointMod$Point
p<-new( Point, 2, 1 )

p$X
p$Y
p$norm()

