# Lattice class ------------------------------------------------------------------------------------
#' @title cflattice
#' @name cflattice
#' @rdname cflattice
#' @aliases cflattice, cflattice-class
#' @description Lattice employed for multinomial valuation
#' @slot lattice lattice of values
#' @slot order order of the lattice
#' @author Pedro Guarderas
#' \email{pedro.felipe.guarderas@@gmail.com}
#' @importFrom methods .valueClassTest new
#' @export
setClass( 'cflattice', slots = representation( lattice = 'list', order = 'numeric' ) )

# Extraction method --------------------------------------------------------------------------------
#' @name [
#' @aliases [, cflattice-methods
#' @docType methods
#' @rdname cflattice
#' @description Extract method for cflattice
#' @param x object of class cflattice
#' @param i index
#' @param j index
#' @param ... extra parameters
#' @param drop drop clause
#' @return Numeric vector with the respective values
#' @author Pedro Guarderas
#' \email{pedro.felipe.guarderas@@gmail.com}
#' @export
setMethod( 
  "[", 
  signature( x = "cflattice", i = "ANY", j = "ANY" ),
  function( x, i, j, ..., drop ) {
    return( x@lattice[[ i, j, ..., drop ]] )
  }
)

# Length method ------------------------------------------------------------------------------------
#' @name length
#' @aliases length, cflattice-methods
#' @docType methods
#' @rdname cflattice
#' @description Extract the length of the cflattice
#' @param x object of class cflattice
#' @return Numeric value
#' @author Pedro Guarderas
#' \email{pedro.felipe.guarderas@@gmail.com}
#' @export
setMethod( 
  "length", 
  signature( x = "cflattice" ),
  function( x ) {
    return( length( x@lattice ) )
  }
)

# Generic order method -----------------------------------------------------------------------------
#' @name order-generic
#' @aliases order, cflattice-methods
#' @docType methods
#' @rdname cflattice
#' @description Extract the order of cflattice
#' @param x object of class cflattice
#' @return Numeric value
#' @author Pedro Guarderas
#' \email{pedro.felipe.guarderas@@gmail.com}
#' @export
setGeneric( 
  "order", 
  valueClass = "numeric", 
  function( x ) {
    standardGeneric( "order" )
  }
)

# Order method -------------------------------------------------------------------------------------
#' @name order
#' @aliases order, cflattice-methods
#' @docType methods
#' @rdname cflattice
#' @description Extract the order of cflattice
#' @param x object of class cflattice
#' @return Numeric value
#' @author Pedro Guarderas
#' \email{pedro.felipe.guarderas@@gmail.com}
#' @export
setMethod( 
  "order", 
  signature( x = "cflattice" ),
  function( x ) {
    return( x@order )
  }
)
