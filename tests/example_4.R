library( Rcpp )
library( RcppArmadillo )
library( CFINI )

sigma<-1
r<-0.05
alpha<-( sigma^2 - 2 * r ) / ( 2 * sigma^2 )
beta<- -( ( 2 * r + sigma^2 ) / ( 2 * sigma^2 ) )^2

t0<-0
t1<-1
t1<-t1 * ( sigma^2 ) / 2
Nt<-500
t<-GridUniform( t0, t1, Nt )
# t<-GridExpAddapt( -3, t0, t1, Nt, 1 )

x0<--8
x1<-8
Nx<-150
x<-GridExpAddapt( -4, x0, x1, Nx, 1 )

K<-1500
I<-sapply( x, FUN = function( x ) exp( -alpha * x ) * max( exp( x ) - K, 0 ) )
# I<-sapply( x, FUN = function( x ) exp( -alpha * x ) * max( K - exp( x ), 0 ) )
A<-rep( 0, Nt )
B<-sapply( t, FUN = function( t ) exp( -alpha * x1 ) * (  exp( x1 ) - K ) * exp( -beta * t ) )
# B<-sapply( t, FUN = function( t ) exp( -alpha * x1 ) * ( K - eSxp( x1 ) ) * exp( -beta * t ) )

print( Nx^2 * sigma^2 *( t1 - t0 ) / ( ( x1 - x0 )^2 ) )

theta<-0.5
U<-DiffusionSolverCNS( 1, theta, I, A, B, t, x )

u<-matrix( 0, Nt, Nx )
for ( n in 1:Nt ) {
  u[ n, ]<-exp( alpha * x ) * exp( beta * ( t1 - t[Nt-n+1] ) ) * U$u[ Nt - n + 1, ]
}
S<-exp( x )
# X11()
persp( t, S, u, theta = 45, phi = 10, xlab = 't', ylab = 'x', zlab = 'u', col = 'olivedrab1',
       box = TRUE, axes = TRUE, main = 'Diffusion solution' )

plot( S, u[Nt,], type = 'l', col = 'red', lwd = 3 )
for ( n in seq( Nt-1, 1, -5 ) ) {
  points( S, u[n,], type = 'l', col = 'red' )
}

