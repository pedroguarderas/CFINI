
#include "cf_ode_solver.h"

//--------------------------------------------------------------------------------------------------
List cf_edo_solv_precor( const Eigen::VectorXd& t,
                         const Eigen::VectorXd& v0,
                         Function f,
                         const int& m,
                         const double& err ) {
  
  int N, d, n, k;
  double dt, nk;
  
  // Grid size
  N = t.size();
  
  // vector dimension
  d = v0.size();
  
  Eigen::VectorXd vk( d ), vk0( d ), fn( d ), fnk( d );
  Eigen::MatrixXd v( N, d );
  
  // Setting conditions
  v.row( 0 ) = v0;
  vk = v0;
  vk0 = vk;
  
  // Predictor corrector solver
  for ( n = 0; n < N - 1; n++ ) {
    dt = t( n + 1 ) - t( n );
    
    // Predictor-corrector step
    fn = Rcpp::as< Eigen::Map< Eigen::VectorXd > >( f( t( n ), vk ) );
    vk = v.row( 0 ).transpose() + dt * fn;
    nk = ( vk - vk0 ).norm() / vk.norm();
    k = 0;
    while ( k <= m || nk > err ) {
      fn = Rcpp::as< Eigen::Map< Eigen::VectorXd > >( f( t( n ), v.row( n ).transpose() ) );
      fnk = Rcpp::as< Eigen::Map< Eigen::VectorXd > >( f( t( n + 1 ), vk ) );
      vk0 = vk;
      vk = v.row( n ).transpose() + 0.5 * dt * ( fn + fnk );
      nk = ( vk - vk0 ).norm() / vk.norm();
      k++;
    }
    v.row( n + 1 ) = vk;
    
  }
  
  return List::create( Named( "t" ) = t, Named( "v" ) = v );
}
