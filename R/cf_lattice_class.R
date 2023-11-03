# Lattice class ------------------------------------------------------------------------------------
#' @title cflattice class
#' @name cflattice-class
#' @rdname cflattice-class
#' @aliases cflattice
#' @description Lattice employed for multinomial valuation
#' @slot lattice lattice of values
#' @slot order order of the lattice
#' @author Pedro Guarderas
#' \email{pedro.felipe.guarderas@@gmail.com}
#' @importFrom methods .valueClassTest new
#' @export
setClass( 'cflattice', slots = representation( lattice = 'list', order = 'numeric' ) )

# Extraction method --------------------------------------------------------------------------------
#' @title [ extraction operator
#' @name [
#' @aliases [
#' @description Extract method for cflattice
#' @param x object of class cflattice
#' @param i index
#' @return Numeric vector with the respective values
#' @author Pedro Guarderas
#' \email{pedro.felipe.guarderas@@gmail.com}
#' @export
setMethod( 
  "[", 
  "cflattice",
  function( x, i ) {
    return( x@lattice[[ i ]] )
  }
)

# Length method ------------------------------------------------------------------------------------
#' @title length
#' @name length
#' @aliases length
#' @description Extract the length of the cflattice
#' @param x object of class cflattice
#' @return Numeric value
#' @author Pedro Guarderas
#' \email{pedro.felipe.guarderas@@gmail.com}
#' @export
setMethod( 
  "length", 
  "cflattice",
  function( x ) {
    return( length( x@lattice ) )
  }
)

# Generic order method -----------------------------------------------------------------------------
#' @title order generic
#' @name order
#' @aliases order
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
#' @title order
#' @name order
#' @aliases order
#' @description Extract the order of cflattice
#' @param x object of class cflattice
#' @return Numeric value
#' @author Pedro Guarderas
#' \email{pedro.felipe.guarderas@@gmail.com}
#' @export
setMethod( 
  "order", 
  "cflattice",
  function( x ) {
    return( x@order )
  }
)
