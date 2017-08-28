library( CFINI )

t0<-0
t1<-2*pi
Nt<-50
t<-GridUniform( t0, t1, Nt )

V0<-c( 1, 0 )
b<-matrix( 0, Nt, 2 )
B<-array( 0, c( 2, 2, Nt ) )
for( n in 1:Nt ) B[,,n]<-matrix( c( 0, -1, 1, 0 ), 2, 2 )

V<-CFThieleSolv( t, V0, b, B )
plot( V$V[,1], V$V[,2], col = 'red', type = 'o', lwd = 2 )
