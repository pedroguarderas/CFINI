
w<-12
nf<-30 * w
n<-seq( 0, nf, 1 ) / w
N<-length(n)

r<-rep( 1 + 0.03 / w, N )
s<-( 1 + 0.02 / w )^n
P<-250 * s
q<-pmin( 1, 0.0005 * exp( seq( 0, 7.7, length.out = N ) ) )
p<-1-q
u<--log(p)
uu<--u

b<-0
for ( i in 1:(N-1) ) {
  b<-c( b, b[i] * r[i] + P[i] * r[i] )
}

R<-1/( 1 + 0.03 / w )
v<-R^(-n)
V<-v * rev( cumsum( P ) ) * rev( cumprod( p ) ) -  R * v * rev( b ) * ( 1 - rev( cumprod( p ) ) )

V1<-0
for ( i in 1:(N-1) ) {
  V1<-c( V1, ( V1[i] + P[i] - R * b[i+1] * q[i] ) / ( R * p[i] ) )
}

plot( n, V, type = 'l', lwd = 2, col = 'midnightblue', panel.first = {abline( h = 0, col = 'red', lwd = 2 )} )
points( n, V1, type = 'l', lwd = 2, col = 'olivedrab3' )



