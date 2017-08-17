#include <Rcpp.h>
using namespace Rcpp;

class Point {
public:
  Point( double x, double y );
  
  double norm() const;
  
  double X, Y;
};

RCPP_MODULE( PointMod ) {
  class_< Point >( "Point" )
  .constructor<double,double>()
  .field( "X", &Point::X )
  .field( "Y", &Point::Y )
  .method( "norm", &Point::norm )
  ;
}
// // You can include R code blocks in C++ files processed with sourceCpp
// // (useful for testing and development). The R code will be automatically 
// // run after the compilation.
// //
// /*** R
// timesTwo(42)
// */
