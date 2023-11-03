setClass( 'cflattice', representation( lattice = 'list' ) )

setMethod( "[", "cflattice",
           function( x, i ) {
             return( x@lattice[[ i ]] )
           }
)
