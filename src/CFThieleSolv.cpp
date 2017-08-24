
#include <RcppArmadillo.h>
#include <map>

using namespace Rcpp;


//' @title Thiele equations solver
//' @description Solver for Thiele equation and computation of Mathematical reserves. The solver 
//' is implemented with a Crank-Nicolson algorithm.
//' @param x0 intial homegity class
//' @param t time grid, could be adapted
//' @param x homegenity classes
//' @param V0 initial mathematical reserve
//' @param r vector with interest rates
//' @param b fixed benefits
//' @param B transition benefits
//' @param u list of inmediate transition matrix
//' @return Return a list with the mathematical reserve, the time grid.
//' @note The solver is implemented to compute the solution forwards in time.
//' @author Pedro Guarderas
//' @export
// [[Rcpp::export]]
List CFThieleSolv( const int& x0,
                   const arma::colvec& t,
                   const arma::colvec& x,
                   const arma::colvec& V0,
                   const arma::colvec& r,
                   const arma::colvec& b,
                   const arma::mat& B,
                   const List& u,
                   const double& theta = 0.5 ) {
  
  int M, N, n, i;
  double dt;
  M = t.size();
  N = V0.size();
  
  arma::mat V( M, N );

  // Setting final conditions  
  for ( i = 0; i < N; i++ ) {
    V( M - 1, i ) = V0( i );
  }
  
  // Thiele equation solver
  for ( n = ( M - 2); n >= 0; n-- ) {
    for ( i = 0; i < N; i++ ){
      dt = t( n + 1 ) - t( n );
      
      V( n, i )<-V( n + 1, i ) - 
          dt * ( r( n + 1 ) * V( n + 1, i ) - b(i) - 
          sum( B[[i]]$b * U[[i]]$u ) - sum( v[B[[i]]$e] * U[[i]]$u ) + 
          V( n + 1, i ) * sum( U[[i]]$u ) );
      }
    }
  }
  
  return  List::create( Named( "V" ) = V,
                        Named( "t" ) = t );
}
    