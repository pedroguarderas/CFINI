
#include "cf_thiele_solv.h"

//--------------------------------------------------------------------------------------------------
// List cf_thiele_solv( const Eigen::VectorXd& t,
//                      const Eigen::VectorXd& V0,
//                      const Eigen::MatrixXd& b,
//                      const Eigen::MatrixXd& B,
//                      const double& theta ) {
// 
//   int M, N, n, i;
//   double dt;
// 
//   N = t.size();
//   M = V0.size();
// 
//   Eigen::VectorXd U( M );
//   Eigen::MatrixXd V( N, M );
//   Eigen::MatrixXd A( M, M );
// 
//   // Setting final conditions
//   V.row( 0 ) = V0.t();
// 
//   // Thiele equation solver
//   for ( n = 0; n < N - 1; n++ ) {
//     dt = t( n + 1 ) - t( n );
//     A.setIdentity();
//     A = A - ( 1.0 - theta ) * dt * B.slice( n + 1 );
//     U = V.row( n ) + theta * dt * ( b.row( n ) + V.row( n ) * B.slice( n ) ) + ( 1.0 - theta ) * dt * b.row( n + 1 );
//     U = arma::trans( arma::solve( A.t(), U.t() ) );
//     V.row( n + 1 ) = U;
//   }
// 
//   return  List::create( Named( "V" ) = V,
//                         Named( "t" ) = t );
// }
//     