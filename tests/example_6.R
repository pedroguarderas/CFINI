

nf<-30 
n<-seq( 0, nf, 1 )
N<-length( n )

r<-rep( 1 + 0.033, N )
s<-( 1 + 0.02 )^n
P<-12 * 250 * s
q<-pmin( 1, 0.0005 * exp( seq( 0, 7.6, length.out = N ) ) )
p<-1-q
u<--log(p)
uu<--u

b<-0
for ( i in 1:(N-1) ) {
  b<-c( b, b[i] * r[i] + P[i] * r[i] )
}

I<-( 1 + 0.033 )
R<-1/I

V<-0
for ( i in 1:(N-1) ) {
  W<-( ( V[i] + P[i] ) * I - b[i+1] * q[i] ) / p[i]
  V<-c( V, W )
}

plot( n, V, type = 'l', lwd = 2, col = 'olivedrab3', ylim = c( 0, 2e5 ),
      panel.first = {abline( h = 0, col = 'red', lwd = 2 )} )

