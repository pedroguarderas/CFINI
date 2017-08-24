library( CFINI )
library( plotly )
library( Rcpp )
library( RcppArmadillo )

sigma<-0.7
r<-0.05
alpha<-( sigma^2 - 2 * r ) / ( 2 * sigma^2 )
beta<- -( ( 2 * r + sigma^2 ) / ( 2 * sigma^2 ) )^2

t0<-0
t1<-1
t1<-t1 * ( sigma^2 ) / 2
Nt<-300
t<-GridUniform( t0, t1, Nt )
# t<-GridExpAddapt( -3, t0, t1, Nt, 1 )

x0<--8
x1<-8
Nx<-100
x<-GridExpAddapt( -4, x0, x1, Nx, 1 )

K<-1500
I<-sapply( x, FUN = function( x ) exp( -alpha * x ) * max( exp( x ) - K, 0 ) )
# I<-sapply( x, FUN = function( x ) exp( -alpha * x ) * max( K - exp( x ), 0 ) )
A<-rep( 0, Nt )
B<-sapply( t, FUN = function( t ) exp( -alpha * x1 ) * (  exp( x1 ) - K ) * exp( -beta * t ) )
# B<-sapply( t, FUN = function( t ) exp( -alpha * x1 ) * ( K - eSxp( x1 ) ) * exp( -beta * t ) )

# print( Nx^2 * sigma^2 *( t1 - t0 ) / ( ( x1 - x0 )^2 ) )

theta<-0.5
# U<-CFDiffSolvCNS( 1, theta, I, A, B, t, x )
V<-CFBlackScholesSolvCNS( sigma, r, theta, I, A, B, t, x )

# u<-matrix( 0, Nt, Nx )
# for ( n in 1:Nt ) {
#   u[ n, ]<-exp( alpha * x ) * exp( beta * ( t1 - t[Nt-n+1] ) ) * U$u[ Nt - n + 1, ]
# }

# plot_ly( x = V$t[,1], y = V$x[,1], z = V$u, alpha = 0.8 ) %>% add_surface()

# View( V$u )
persp( x = V$t[,1], y = V$x[,1], z = V$u, col = 'red',
       theta = 40, xlab = 't', ylab = 'S', zlab = 'V' )
# persp( t, S, u )

plot( V$x[,1], V$u[Nt,], type = 'l', col = 'red',
      xlab = 'S', ylab = 'u',
      xlim = c( 0, 3000 ), ylim = c( 0, 1750 ), axes = FALSE )
axis( 1, at=seq( 0, 3000, 250 ), labels=seq( 0, 3000, 250 ) )
axis( 2, at=seq( 0, 1750, 250 ), labels=seq( 0, 1750, 250 ) )
for ( n in seq( Nt-1, 1, -5 ) ) {
  points( V$x[,1], V$u[n,], type = 'l', col = 'red',
        xlab = 'S', ylab = 'u',
        xlim = c( 0, 3000 ), ylim = c( 0, 1750 ) )
}

# library( animation )
# saveVideo( {
#   plot( S, u[Nt,], type = 'l', col = 'red', 
#         xlab = 'S', ylab = 'u', 
#         xlim = c( 0, 3000 ), lim = c( 0, 1750 ), axes = FALSE )
#   axis( 1, at=seq( 0, 3000, 250 ), labels=seq( 0, 3000, 250 ) )
#   axis( 2, at=seq( 0, 1750, 250 ), labels=seq( 0, 1750, 250 ) )
#   for ( n in seq( Nt-1, 1, -1 ) ) {
#     plot( S, u[n,], type = 'l', col = 'red', 
#           xlab = 'S', ylab = 'u',
#           xlim = c( 0, 3000 ), ylim = c( 0, 1750 ) )
#     axis( 1, at=seq( 0, 3000, 250 ), labels=seq( 0, 3000, 250 ) )
#     axis( 2, at=seq( 0, 1750, 250 ), labels=seq( 0, 1750, 250 ) )
#   }
# }, ffmpeg = '/usr/bin/ffmpeg', 
# interval = 0.01, video.name = 'BlackScholesAni.mp4', ani.width = 600, ani.height = 600 )

  