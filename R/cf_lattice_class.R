# Lattice class ------------------------------------------------------------------------------------
#' An S4 class for lattice
#' @rdname cflattice_class
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
#' Extraction method for cflattice
#' @rdname cflattice_class
#' @aliases [, ANY, ANY-method
#' @docType methods
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
#' Get length of cflattice
#' @rdname cflattice_class
#' @aliases length, numeric, ANY-method
#' @docType methods
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
#' Generic method for order of the cflattice
#' @rdname cflattice_class
#' @aliases order, numeric, ANY-method
#' @docType methods
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
#' Method for order of the cflattice
#' @rdname cflattice_class
#' @aliases order, numeric, ANY-method
#' @docType methods
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
