
#include "cf_diff_solv.h"

//--------------------------------------------------------------------------------------------------
List cf_diff_solv_euls( const double& alpha,
                        const Eigen::VectorXd& I,
                        const Eigen::VectorXd& A,
                        const Eigen::VectorXd& B,
                        const Eigen::VectorXd& t,
                        const Eigen::VectorXd& x ) {
  
  int n, i;
  double dt, dxf, dxb, lambda, h;
  int Nt, Nx;
  
  Nt = t.size();
  Nx = x.size();
  Eigen::MatrixXd u( Nt, Nx );
  
  // Inicial condition
  for ( i = 0; i < Nx; i++ ) {
    u( 0, i ) = I( i );
  }
  
  // Boundary conditions
  for ( n = 0; n < Nt; n++ ) {
    u( n, 0 ) = A( n );
    u( n, Nt - 1 ) = B( n );
  }
  
  // Solver
  // Euler explicit scheme
  for ( n = 0; n < Nt - 1; n++ ) {
    dt = t( n + 1 ) - t( n );
    for ( i = 1; i < Nx - 1; i++ ) {
      dxf = x( i + 1 ) - x( i );
      dxb = x( i ) - x( i - 1 );
      lambda = alpha * dt / ( dxf * dxb );
      u( n + 1, i ) = u( n, i ) + lambda * ( u( n, i + 1 ) - 2.0 * u( n, i ) + u( n, i - 1 ) );
    }
  }
  
  return List::create( Named( "u" ) = u,
                       Named( "t" ) = t,
                       Named( "x" ) = x );
}

//--------------------------------------------------------------------------------------------------
List cf_diff_solv_cns( const double& alpha,
                       const double& theta,
                       const Eigen::VectorXd& I,
                       const Eigen::VectorXd& A,
                       const Eigen::VectorXd& B,
                       const Eigen::VectorXd& t,
                       const Eigen::VectorXd& x ) {
  
  int n, i;
  double dt, dxf, dxb, lambdaf, lambdab, h;
  int Nt, Nx, nt, nx;
  Nt = t.size();
  Nx = x.size();
  
  nx = Nx - 1;
  nt = Nt - 1;
  
  Eigen::VectorXd a( Nx );
  Eigen::VectorXd b( Nx );
  Eigen::VectorXd c( Nx );
  Eigen::VectorXd d( Nx );
  Eigen::MatrixXd u( Nt, Nx );
  
  // Setting the inicial condition
  for ( i = 0; i < Nx; i++ ) {
    u( 0, i ) = I( i );
  }
  
  // Setting boundary conditions
  for ( n = 0; n < Nt; n++ ) {
    u( n, 0 ) = A( n );
    u( n, nx ) = B( n );
  }
  
  // Solver
  // Crank-Nicolson scheme
  for ( n = 0; n < nt; n++ ) {
    
    // Time difference
    dt = t( n + 1 ) - t( n );
    
    d( 0 ) = A( n );
    d( nx ) = B( n );
    
    // forward difference
    dxf = x( 1 ) - x( 0 );
    lambdaf = alpha * dt / ( dxf * dxf );
    
    // backward difference
    dxb = x( nx ) - x( nx - 1 );
    lambdab = alpha * dt / ( dxb * dxf );
    
    // Initial filling
    a( 0 ) = -theta * lambdab;
    a( nx ) = -theta * lambdab;
    b( 0 ) = 1.0 + theta * ( lambdaf + lambdab );
    b( nx ) = 1.0 + theta * ( lambdaf + lambdab );
    c( 0 ) = -theta * lambdaf;
    c( nx ) = -theta * lambdaf;
    
    for ( i = 1; i < nx; i++ ) {
      dxf = x( i + 1 ) - x( i );
      dxb = x( i ) - x( i - 1 );
      lambdab = alpha * dt / ( dxb * dxf );
      lambdaf = alpha * dt / ( dxf * dxf );
      
      a( i ) = -theta * lambdab;
      b( i ) = 1.0 + theta * ( lambdaf + lambdab );
      c( i ) = -theta * lambdaf;
      
      d( i ) = u( n, i ) + ( 1.0 - theta ) * 
        ( lambdaf * ( u( n, i + 1 ) - u( n, i ) ) - lambdab * ( u( n, i ) - u( n, i - 1 ) ) ) ;
    }
    
    cf_tri_diag_solv( a, b, c, d );
    d( 0 ) = A( n );
    d( nx ) = B( n );
    
    for ( i = 0; i < Nx; i++ ) {
      u( n + 1, i ) = d( i );
    }
    
  }
  
  return List::create( Named( "u" ) = u,
                       Named( "t" ) = t,
                       Named( "x" ) = x );
}

//--------------------------------------------------------------------------------------------------
List cf_black_scholes_solv_cns( const double& sigma,
                                const double& rate,
                                const double& theta,
                                const Eigen::VectorXd& I,
                                const Eigen::VectorXd& A,
                                const Eigen::VectorXd& B,
                                const Eigen::VectorXd& t,
                                const Eigen::VectorXd& x ) {
  
  int n, i;
  double dt, dxf, dxb, lambdaf, lambdab, h;
  double D, alpha, beta;
  int Nt, Nx, nt, nx;
  Nt = t.size();
  Nx = x.size();
  
  D = 1.0;
  alpha = ( sigma * sigma - 2.0 * rate ) / ( 2.0 * sigma * sigma );
  beta = ( 2.0 * rate + sigma * sigma ) / ( 2.0 * sigma * sigma );
  beta = -beta * beta;
  
  nx = Nx - 1;
  nt = Nt - 1;
  
  Eigen::VectorXd a( Nx );
  Eigen::VectorXd b( Nx );
  Eigen::VectorXd c( Nx );
  Eigen::VectorXd d( Nx );
  
  Eigen::MatrixXd u( Nt, Nx );
  Eigen::MatrixXd v( Nt, Nx );
  
  // Inicial condition
  for ( i = 0; i < Nx; i++ ) {
    u( 0, i ) = I( i );
  }
  
  // Boundary conditions
  for ( n = 0; n < Nt; n++ ) {
    u( n, 0 ) = A( n );
    u( n, nx ) = B( n );
  }
  
  // Solver
  // Crank-Nicolson scheme
  for ( n = 0; n < nt; n++ ) {
    dt = t( n + 1 ) - t( n );
    
    d( 0 ) = A( n );
    d( nx ) = B( n );
    
    dxf = x( 1 ) - x( 0 );
    lambdaf = D * dt / ( dxf * dxf );
    
    dxb = x( nx ) - x( nx - 1 );
    lambdab = D * dt / ( dxb * dxf );
    
    a( 0 ) = -theta * lambdab;
    a( nx ) = -theta * lambdab;
    b( 0 ) = 1.0 + theta * ( lambdaf + lambdab );
    b( nx ) = 1.0 + theta * ( lambdaf + lambdab );
    c( 0 ) = -theta * lambdaf;
    c( nx ) = -theta * lambdaf;
    
    for ( i = 1; i < nx; i++ ) {
      dxf = x( i + 1 ) - x( i );
      dxb = x( i ) - x( i - 1 );
      lambdab = D * dt / ( dxb * dxf );
      lambdaf = D * dt / ( dxf * dxf );
      
      a( i ) = -theta * lambdab;
      b( i ) = 1.0 + theta * ( lambdaf + lambdab );
      c( i ) = -theta * lambdaf;
      
      d( i ) = u( n, i ) + ( 1.0 - theta ) * 
        ( lambdaf * ( u( n, i + 1 ) - u( n, i ) ) - lambdab * ( u( n, i ) - u( n, i - 1 ) ) ) ;
    }
    
    cf_tri_diag_solv( a, b, c, d );
    d( 0 ) = A( n );
    d( nx ) = B( n );
    
    for ( i = 0; i < Nx; i++ ) {
      u( n + 1, i ) = d( i );
    }
    
  }
  
  for ( i = 0; i < Nx; i++ ) {
    for ( n = 0; n < Nt; n++ ) {
      v( n, i ) = exp( alpha * x( i ) ) * exp( beta * ( t( nt ) - t( nt - n ) ) ) * u( nt - n, i );
    }
    d( i ) = exp( x(i) );
  }
  
  return List::create( Named( "u" ) = v,
                       Named( "t" ) = t,
                       Named( "x" ) = d );
}
