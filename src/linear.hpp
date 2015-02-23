//__________________________________________________________________________________________________

#include <cmath>
#include <iostream>
#include <vector>

using namespace std;


//__________________________________________________________________________________________________
/* Vector class definition */
template< typename Field = double >
class Vector : public vector< Field > {
public:
	Vector() {}

	Vector( const size_t n ) :
		vector< Field >( n ) {}

	Vector( const vector< Field > v ) :
		vector< Field >( v ) {}

	Vector( Field* begin, Field* end ) :
		vector< Field >( begin, end ) {}

	~Vector() {}

	Vector< Field > add( const size_t i, const Field& x ) const {
		Vector< Field > z( this->size() );

		#pragma omp parallel for

		for( size_t j = 0; j < this->size(); j++ ) {
			z[j] = this->at( j );
		}

		z[i] += x;
		return z;
	}

	friend std::ostream& operator<< ( std::ostream& out, Vector< Field >& x ) {
		for( size_t i = 0; i < x.size(); i++ ) {
			out << x[i] << " ";
		}

		return out;
	}

};

template< typename Field = double >
Field norm( const Vector< Field >& x ) {
	size_t i;
	Field n;

	n = Field( 0.0 );
	#pragma omp parallel for private( i ) reduction (+:n)

	for( i = 0; i < x.size(); i++ ) {
		n = n + x[i] * x[i];
	}

	return sqrt( n );
}

template< typename Field = double >
Vector< Field > operator+ ( const Vector< Field >& x, const Vector< Field >& y ) {
	size_t m, M, i;
	m = min( x.size(), y.size() );
	M = max( x.size(), y.size() );

	Vector< Field > z( M );

	#pragma omp parallel for private(i)

	for( i = 0; i < m; i++ ) {
		z[i] = x[i] + y[i];
	}

	if( x.size() > y.size() ) {
		#pragma omp parallel for private(i)

		for( i = m; i < M; i++ ) {
			z[i] = y[i];
		}
	}
	else if( x.size() < y.size() ) {
		#pragma omp parallel for private(i)

		for( i = m; i < M; i++ ) {
			z[i] = x[i];
		}
	}

	return z;
}

template< typename Field = double >
Vector< Field > operator- ( const Vector< Field >& x, const Vector< Field >& y ) {
	size_t m, M, i;
	m = min( x.size(), y.size() );
	M = max( x.size(), y.size() );

	Vector< Field > z( M );

	#pragma omp parallel for private(i)

	for( i = 0; i < m; i++ ) {
		z[i] = x[i] - y[i];
	}

	if( x.size() > y.size() ) {
		#pragma omp parallel for private(i)

		for( i = m; i < M; i++ ) {
			z[i] = y[i];
		}
	}
	else if( x.size() < y.size() ) {
		#pragma omp parallel for private(i)

		for( i = m; i < M; i++ ) {
			z[i] = x[i];
		}
	}

	return z;
}

template< typename Field = double >
Vector< Field > operator+ ( const Field a, const Vector< Field >& x ) {
	size_t i;
	Vector< Field > z( x.size() );

	#pragma omp parallel for private(i)

	for( i = 0; i < x.size(); i++ ) {
		z[i] = a + x[i];
	}

	return z;
}

template< typename Field = double >
Vector< Field > operator+ ( const Vector< Field >& x, const Field a ) {
	Vector< Field > z( x.size() );

	#pragma omp parallel for

	for( size_t i = 0; i < x.size(); i++ ) {
		z[i] = a + x[i];
	}

	return z;
}

template< typename Field = double >
Vector< Field > operator* ( const Vector< Field >& x, const Vector< Field >& y ) {
	size_t s;

	if( x.size() < y.size() ) {
		s = x.size();
	}
	else {
		s = y.size();
	}

	Vector< Field > z( s );

	#pragma omp parallel for

	for( size_t i = 0; i < s; i++ ) {
		z[i] = x[i] * y[i];
	}

	return z;
}

template< typename Field = double >
Vector< Field > operator/ ( const Vector< Field >& x, const Vector< Field >& y ) {
	bool dzero;
	size_t s;

	if( x.size() < y.size() ) {
		s = x.size();
	}
	else {
		s = y.size();
	}

	Vector< Field > z( s );
	dzero = false;

	//    #pragma omp parallel for
	for( size_t i = 0; i < s; i++ ) {
		if( y[i] != Field( 0.0 ) ) {
			z[i] = x[i] / y[i];
		}
		else {
			dzero = true;
			break;
		}
	}

	if( dzero ) {
		z.clear();
	}

	return z;
}

template< typename Field = double >
Vector< Field > operator* ( const Field a, const Vector< Field >& x ) {
	Vector< Field > z( x.size() );

	#pragma omp parallel for

	for( size_t i = 0; i < x.size(); i++ ) {
		z[i] = a * x[i];
	}

	return z;
}

template< typename Field = double >
Vector< Field > operator* ( const Vector< Field >& x, const Field a ) {
	Vector< Field > z( x.size() );

	#pragma omp parallel for

	for( size_t i = 0; i < x.size(); i++ ) {
		z[i] = a * x[i];
	}

	return z;
}

template< typename Field = double >
Field scalar( const Vector< Field >& x, const Vector< Field >& y ) {
	size_t m, i;
	Field s;

	m = min( x.size(), y.size() );
	s = Field( 0.0 );

	#pragma omp parallel for reduction (+:s)

	for( i = 0; i < m; i++ ) {
		s = s + x[i] * y[i];
	}

	return s;
}


//__________________________________________________________________________________________________
/* Matrix class definition */
template< typename Field = double >
class Matrix : public vector< Field > {
public:
	Matrix() {}

	Matrix( const size_t n ) :
		vector< Field >( n ) {
	}

	Matrix( const size_t m, const size_t n ) :
		vector< Field >( m* n ),
		_dim( { m, n } ) {
	}

	~Matrix() {}

	Field get( size_t i, size_t j ) const {
		return this->at( j + i * _dim[1] );
	}

	size_t dim( size_t i ) const {
		return _dim[i];
	}

	void set_dim( vector< size_t > dim ) {
		_dim = dim;
	}

	// supposed that the matrix has already elements
	void set_column( const size_t j, const Vector< Field >& x ) {
		if( x.size() == _dim[0] && j < _dim[1] ) {
			#pragma omp parallel for

			for( size_t i = 0; i < _dim[0]; i++ ) {
				this->at( j + i * _dim[1] ) = x[i];
			}
		}
	}

	// supposed that the matrix has already elements
	void set_row( const size_t i, const Vector< Field >& x ) {
		if( x.size() == _dim[1]  && i < _dim[0] ) {
			#pragma omp parallel for

			for( size_t j = 0; j < _dim[1]; j++ ) {
				this->at( j + i * _dim[1] ) = x[j];
			}
		}
	}

	friend std::ostream& operator<< ( std::ostream& out, Matrix< Field >& x ) {
		for( size_t i = 0; i < x.dim( 0 ); i++ ) {
			for( size_t j = 0; j < x.dim( 1 ); j++ ) {
				out << x.get( i, j ) << " ";
			}

			if( i < x.dim( 0 ) - 1 ) {
				out << endl;
			}
		}

		return out;
	}
	
	friend Field operator[]( const size_t i, const size_t j ) {
    return this->at( j + i * _dim[1] );
  }

private:
	vector< size_t > _dim;
};

template< typename Field = double >
Matrix< Field > operator+ ( const Matrix< Field >& x, const Matrix< Field >& y ) {
	Matrix< Field > z;
	z.set_dim( { x.dim( 0 ), x.dim( 1 ) } );

	if( x.size() == y.size() ) {
		z.resize( x.size() );

		#pragma omp parallel shared(x,y,z) for

		for( size_t i = 0; i < x.size(); i++ ) {
			z[i] = x[i] + y[i];
		}
	}

	return z;
}

template< typename Field = double >
Matrix< Field > operator- ( const Matrix< Field >& x, const Matrix< Field >& y ) {
	Matrix< Field > z;
	z.set_dim( { x.dim( 0 ), x.dim( 1 ) } );

	if( x.size() == y.size() ) {
		z.resize( x.size() );

		#pragma omp parallel shared(x,y,z) for

		for( size_t i = 0; i < x.size(); i++ ) {
			z[i] = x[i] - y[i];
		}
	}

	return z;
}

template< typename Field = double >
Matrix< Field > operator* ( Matrix< Field >& x, Matrix< Field >& y ) {
	Field sum ;
	Matrix< Field > z;

	if( x.dim( 1 ) == y.dim( 0 ) ) {

		z.resize( x.dim( 0 ) * y.dim( 1 ) );
		z.set_dim( { x.dim( 0 ), y.dim( 1 ) } );

		#pragma omp parallel shared(x,y,z) for collapse(2)

		for( size_t i = 0; i < x.dim( 0 ); i++ ) {
			for( size_t k = 0; k < y.dim( 1 ); k++ ) {
				for( size_t j = 0; j < x.dim( 1 ); j++ ) {
					if( j == 0 ) {
						z[i, k] = 0.0;
					}

					z[i, k] = z[i, k] + x.get( i, j ) * y.get( j, k );
				}

				z.push_back( sum );
			}
		}
	}

	return z;
}

template< typename Field = double >
Vector< Field > operator* ( Matrix< Field >& x, Vector< Field >& y ) {
	Field sum;
	Vector< Field > z( x.dim( 0 ) );

	if( x.dim( 1 ) == y.size() ) {

		#pragma omp parallel for

		for( size_t i = 0; i < x.dim( 0 ); i++ ) {
			z[i] = Field( 0.0 );
		}

		#pragma omp parallel for collapse(2)

		for( size_t i = 0; i < x.dim( 0 ); i++ ) {
			for( size_t j = 0; j < x.dim( 1 ); j++ ) {
				z[i] = z[i] + x.get( i, j ) * y[j];
			}
		}
	}

	return z;
}

template< typename Field = double >
Vector< Field > operator* ( Vector< Field >& y, Matrix< Field >& x ) {
	Field sum;
	Vector< Field > z( x.dim( 1 ) );

	if( x.dim( 0 ) == y.size() ) {

		#pragma omp parallel for

		for( size_t j = 0; j < x.dim( 1 ); j++ ) {
			z[j] = Field( 0.0 );
		}

		#pragma omp parallel for collapse(2)

		for( size_t i = 0; i < x.dim( 1 ); i++ ) {
			for( size_t j = 0; j < x.dim( 0 ); j++ ) {
				z[i] = z[i] + x.get( j, i ) * y[j];
			}
		}
	}

	return z;
}

template< typename Field = double >
Matrix< Field > operator+ ( const Field a, const Matrix< Field >& x ) {
	Matrix< Field > z( x.size() );
	z.set_dim( { x.dim( 0 ), x.dim( 1 ) } );

	#pragma omp parallel for

	for( size_t i = 0; i < x.size(); i++ ) {
		z[i] = a + x[i];
	}

	return z;
}

template< typename Field = double >
Matrix< Field > operator+ ( const Matrix< Field >& x, const Field a ) {
	Matrix< Field > z( x.size() );
	z.set_dim( { x.dim( 0 ), x.dim( 1 ) } );

	#pragma omp parallel for

	for( size_t i = 0; i < x.size(); i++ ) {
		z[i] = a + x[i];
	}

	return z;
}

template< typename Field = double >
Matrix< Field > operator* ( const Field a, const Matrix< Field >& x ) {
	Matrix< Field > z( x.size() );
	z.set_dim( { x.dim( 0 ), x.dim( 1 ) } );

	#pragma omp parallel for

	for( size_t i = 0; i < x.size(); i++ ) {
		z[i] = a * x[i];
	}

	return z;
}

template< typename Field = double >
Matrix< Field > operator* ( const Matrix< Field >& x, const Field a ) {
	Matrix< Field > z( x.size() );
	z.set_dim( { x.dim( 0 ), x.dim( 1 ) } );

	#pragma omp parallel for

	for( size_t i = 0; i < x.size(); i++ ) {
		z[i] = a * x[i];
	}

	return z;
}

template< typename Field = double >
Matrix< Field > Trp( const Matrix< Field >& x ) {
	Matrix< Field > z( x.dim( 0 ) * x.dim( 1 ) );

	z.set_dim( { x.dim( 1 ), x.dim( 0 ) } );

	#pragma omp parallel for

	for( size_t j = 0; j < x.dim( 1 ); j++ ) {
		#pragma omp parallel for

		for( size_t i = 0; i < x.dim( 0 ); i++ ) {
			z[ i + j * x.dim( 0 ) ] = x.get( i, j );
		}
	}

	return z;
}

template< typename Field = double >
Field Tr( const Matrix< Field >& x ) {
	size_t i;
	double m;
	Field tr;

	tr = Field( 0.0 );
	m = x.dim( 0 ) * x.dim( 1 );

	#pragma omp parallel for reduction (+:tr)

	for( i = 0; i < m; i++ ) {
		tr = tr + x.get( i, i );
	}

	return tr;
}

template< typename Field = double >
Field scalar( const Matrix< Field >& x, const Matrix< Field >& y ) {
	Field n;

	n = Field( 0.0 );

	if( x.size() == y.size() ) {
		#pragma omp parallel for reduction (+:n)

		for( size_t i = 0; i < x.size(); i++ ) {
			n = n + x[i] * y[i];
		}
	}

	return n;
}
