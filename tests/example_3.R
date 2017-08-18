library( Rcpp )
library( RcppArmadillo )
library( CFINI )

t0<-0
t1<-0.3
# t1<-0.543
Nt<-100
t<-GridUniform( t0, t1, Nt )

x0<-0
x1<-1
Nx<-100
x<-GridUniform( x0, x1, Nx )

# I<-sapply( x, FUN = function( x ) sin( pi * x ) )
# I<-sapply( x, FUN = function( x ) exp( -( x - 0.5 )^2 ) / sqrt( 2 * pi ) )
I<-sapply( x, FUN = function( x ) if ( x >= 0.4 & x <= 0.6 ) return( 1 ) else return( 0 ) )
A<-rep( 0, Nt )
B<-rep( 0, Nt )

alpha<-1e-2
theta<-0.5

hx<-(x1-x0)/(Nx-1)
ht<-(t1-t0)/(Nt-1)
ht / ( hx * hx )
1/alpha

U<-DiffusionSolverCNS( alpha, theta, I, A, B, t, x )

persp( t, x, U$u, theta = 120, phi = 20, xlab = 't', ylab = 'x', zlab = 'u', col = 'gold',
       box = TRUE, axes = TRUE, main = 'Diffusion solution' )

