
#include "cf_ode_solv.h"

//--------------------------------------------------------------------------------------------------
List cf_edo_solv_precor( const Eigen::VectorXd& t,
                         const Eigen::VectorXd& v0,
                         const Function& f,
                         const int& m,
                         const double& err ) {
  
  int N, d, n, k;
  double dt, nk;
  
  // Grid size
  N = t.size();
  
  // vector dimension
  d = v0.size();
  
  Eigen::VectorXd vk( d ), vk0( d );
  Eigen::MatrixXd v( N, d );
  
  // Setting conditions
  v.row( 0 ) = v0.transpose();
  
  // Predictor corrector solver
  for ( n = 0; n < N - 1; n++ ) {
    dt = t( n + 1 ) - t( n );
    
    // Predictor-corrector step
    vk0 = v.row( n ) + dt * Rcpp::as< Eigen::Map< Eigen::VectorXd > >( f( t( n ), vk ) );
    k = 0;
    while ( k <= m || nk > err ) {
      vk = v.row( n ) 
      + 0.5 * dt * ( Rcpp::as< Eigen::Map< Eigen::VectorXd > >( f( t( n ), v.row( n + 1 ) ) ) 
      + Rcpp::as< Eigen::Map< Eigen::VectorXd > >( f( t( n + 1 ), vk0 ) ) );
      vk0 = vk;
      nk = ( vk - vk0 ).norm() / vk.norm();
      k++;
    }
    v.row( n + 1 ) = vk;
    
  }
  
  return  List::create( Named( "t" ) = t, Named( "v" ) = v );
}
