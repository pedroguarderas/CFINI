#include "point.h"

Point::Point( double x, double y ) : X(x), Y(y) {
}
  
double Point::norm() const {
    return( X * X + Y * Y );
}
