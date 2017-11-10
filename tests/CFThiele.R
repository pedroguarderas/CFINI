library( CFINI )
library( Matrix )

x<-0:105
n<-length(x)
p12<-as.vector( GridExpAddapt( -8, 0.1, 1, n, 1 ) )
p13<-as.vector( GridExpAddapt( 15, 0.999, 0.00001, n, 1 ) )
p23<-as.vector( GridExpAddapt( 5, 0.9999, 0.000001, n, 1 ) )
u12<--log(p12)
u13<--log(p13)
u23<--log(p23)

Q<-list()
for ( i in 1:n ) {
  Q[[i]]<-matrix( c( -u12[i] - u13[i], u12[i], u13[i],
                     0, -u23[i], u23[i],
                     0, 0, 0 ), 3, 3, byrow = TRUE )
}

P<-list()
for ( i in 1:n ) {
  E<-eigen( Q[[i]] )
  D<-diag( exp( E$values ) )
  P[[i]]<-E$vectors %*% D %*% solve( E$vectors )
}

y<-33
s<-1
t<-40
Pr<-list( P[[y+s]] )
i<-1
for ( k in s:t ) {
  Pr[[i+1]]<-Pr[[i]] %*% P[[y+k]]
  i<-i+1
}

P[[33]]
Pr[[10]]

sapply( P, FUN = rowSums )
apply( sapply( Q, FUN = rowSums ), c(1,2), FUN = function(x) ifelse( x < 1e-15, 0, x ) )

plot( 1-p13 )
points( 1-p12 )

plot( 1-p23)
