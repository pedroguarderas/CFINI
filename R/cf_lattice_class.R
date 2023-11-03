setClass( 'cflattice', representation( lattice = 'list', order = 'numeric' ) )

setMethod( 
  "[", 
  "cflattice",
  function( x, i ) {
    return( x@lattice[[ i ]] )
  }
)

setMethod( 
  "length", 
  "cflattice",
  function( x ) {
    return( length( x@lattice ) )
  }
)

setGeneric( 
  "order", 
  valueClass = "numeric", 
  function( x ) {
    standardGeneric( "order" )
  }
)

setMethod( 
  "order", 
  "cflattice",
  function( x ) {
    return( x@order )
  }
)
