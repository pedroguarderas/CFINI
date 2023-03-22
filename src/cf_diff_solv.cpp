
#include "cf_diff_solv.h"

//--------------------------------------------------------------------------------------------------
List cf_diff_solv_euls( const Eigen::MatrixXd& alpha,
                        const Eigen::VectorXd& u0,
                        const Eigen::VectorXd& u1,
                        const Eigen::VectorXd& u2,
                        const Eigen::VectorXd& t,
                        const Eigen::VectorXd& x,
                        const bool is_initial = true ) {
  
  int n, i;
  double dt, dxf, dxb, lambda;
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

  // Setting boundary conditions
  u.col( 0 ) = u1;
  u.col( nx ) = u2;
  
  // Setting the inicial condition
  if ( is_initial ) {
    u.row( 0 ) = u0;
    
    // Solver
    // Euler implicit scheme
    for ( n = 0; n < nt; n++ ) {
      
      // Time difference
      dt = t( n + 1 ) - t( n );
      
      d( 0 ) = u1( n );
      d( nx ) = u2( n );
      
      // Space forward difference
      dxf = x( 1 ) - x( 0 );
      lambda = alpha( n + 1, 0 ) * dt / ( dxf * dxf );
      
      // Initial filling
      a( 0 ) = -lambda;
      a( nx ) = -lambda;
      b( 0 ) = 1.0 + 2.0 * lambda;
      b( nx ) = 1.0 + 2.0 * lambda;
      c( 0 ) = -lambda;
      c( nx ) = -lambda;
      
      for ( i = 1; i < nx; i++ ) {
        
        dxf = x( i + 1 ) - x( i );
        dxb = x( i ) - x( i - 1 );
        
        lambda = alpha( n + 1, i ) * dt / ( dxb * dxf );
        
        a( i ) = -lambda;
        b( i ) = 1.0 + 2.0 * lambda;
        c( i ) = -lambda;
        d( i ) = u( n, i );
      }
      
      cf_tri_diag_solv( a, b, c, d );
      d( 0 ) = u1( n );
      d( nx ) = u2( n );
      u.row( n + 1 ) = d;
      
    }
  } else {
    u.row( nt ) = u0;
    
    // Solver
    // Euler implicit scheme
    for ( n = nt - 1; n >= 0; n-- ) {
      
      // Time difference
      dt = t( n + 1 ) - t( n );
      
      d( 0 ) = u1( n );
      d( nx ) = u2( n );
      
      // Space forward difference
      dxf = x( 1 ) - x( 0 );
      lambda = alpha( n, 0 ) * dt / ( dxf * dxf );
      
      // Initial filling
      a( 0 ) = -lambda;
      a( nx ) = -lambda;
      b( 0 ) = 1.0 + 2.0 * lambda;
      b( nx ) = 1.0 + 2.0 * lambda;
      c( 0 ) = -lambda;
      c( nx ) = -lambda;
      
      for ( i = 1; i < nx; i++ ) {
        
        dxf = x( i + 1 ) - x( i );
        dxb = x( i ) - x( i - 1 );
        
        lambda = alpha( n, i ) * dt / ( dxb * dxf );
        
        a( i ) = -lambda;
        b( i ) = 1.0 + 2.0 * lambda;
        c( i ) = -lambda;
        d( i ) = u( n + 1, i );
      }
      
      cf_tri_diag_solv( a, b, c, d );
      d( 0 ) = u1( n );
      d( nx ) = u2( n );
      u.row( n ) = d;
      
    }
  }
  
  return List::create( Named( "u" ) = u,
                       Named( "t" ) = t,
                       Named( "x" ) = x );
}

//--------------------------------------------------------------------------------------------------
List cf_diff_solv_cns( const double& theta,
                       const Eigen::MatrixXd& alpha,
                       const Eigen::VectorXd& u0,
                       const Eigen::VectorXd& u1,
                       const Eigen::VectorXd& u2,
                       const Eigen::VectorXd& t,
                       const Eigen::VectorXd& x,
                       const bool is_initial = true ) {
  
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
  
  // Setting boundary conditions
  u.col( 0 ) = u1;
  u.col( nx ) = u2;
  
  // Setting the inicial condition
  if ( is_initial ) {
    u.row( 0 ) = u0;
    
    // Solver
    // Crank-Nicolson scheme
    for ( n = 0; n < nt; n++ ) {
      
      // Time difference
      dt = t( n + 1 ) - t( n );
      
      d( 0 ) = u1( n );
      d( nx ) = u2( n );
      
      // forward difference
      dxf = x( 1 ) - x( 0 );
      lambdaf = alpha( n + 1, 0 ) * dt / ( dxf * dxf );
      
      // backward difference
      dxb = x( nx ) - x( nx - 1 );
      lambdab = alpha( n + 1, nx ) * dt / ( dxb * dxf );
      
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
        lambdab = alpha( n + 1, i ) * dt / ( dxb * dxf );
        lambdaf = alpha( n + 1, i ) * dt / ( dxf * dxf );
        
        a( i ) = -theta * lambdab;
        b( i ) = 1.0 + theta * ( lambdaf + lambdab );
        c( i ) = -theta * lambdaf;
        
        d( i ) = u( n, i ) + ( 1.0 - theta ) * 
          ( lambdaf * ( u( n, i + 1 ) - u( n, i ) ) - 
            lambdab * ( u( n, i ) - u( n, i - 1 ) ) );
        
      }
      
      cf_tri_diag_solv( a, b, c, d );
      d( 0 ) = u1( n );
      d( nx ) = u2( n );
      u.row( n + 1 ) = d;
      
    }
  
  } else {
    
    u.row( nt ) = u0;
    
    // Solver
    // Crank-Nicolson scheme
    for ( n = nt - 1; n >= 0; n-- ) {
      
      // Time difference
      dt = t( n + 1 ) - t( n );
      
      d( 0 ) = u1( n );
      d( nx ) = u2( n );
      
      // forward difference
      dxf = x( 1 ) - x( 0 );
      lambdaf = alpha( n, 0 ) * dt / ( dxf * dxf );
      
      // backward difference
      dxb = x( nx ) - x( nx - 1 );
      lambdab = alpha( n, nx ) * dt / ( dxb * dxf );
      
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
        lambdab = alpha( n, i ) * dt / ( dxb * dxf );
        lambdaf = alpha( n, i ) * dt / ( dxf * dxf );
        
        a( i ) = -theta * lambdab;
        b( i ) = 1.0 + theta * ( lambdaf + lambdab );
        c( i ) = -theta * lambdaf;
        
        d( i ) = u( n + 1, i ) + ( 1.0 - theta ) * 
          ( lambdaf * ( u( n + 1, i + 1 ) - u( n + 1, i ) ) - 
            lambdab * ( u( n + 1, i ) - u( n + 1, i - 1 ) ) );
      }
      
      cf_tri_diag_solv( a, b, c, d );
      d( 0 ) = u1( n );
      d( nx ) = u2( n );
      u.row( n ) = d;
      
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
                                const Eigen::VectorXd& u0,
                                const Eigen::VectorXd& u1,
                                const Eigen::VectorXd& u2,
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
  
  // Setting boundary conditions
  u.col( 0 ) = u1;
  u.col( nx ) = u2;
  
  // Inicial condition
  for ( i = 0; i < Nx; i++ ) {
    u( 0, i ) = I( i );
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
        ( lambdaf * ( u( n, i + 1 ) - u( n, i ) ) - 
          lambdab * ( u( n, i ) - u( n, i - 1 ) ) ) ;
    }
    
    cf_tri_diag_solv( a, b, c, d );
    d( 0 ) = u1( n );
    d( nx ) = u2( n );
    
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
