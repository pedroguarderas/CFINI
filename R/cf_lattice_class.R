# Lattice class ------------------------------------------------------------------------------------
#' @title Lattice class for multinomial valuation
#' @description Lattice employed for multinomial valuation
#' @return A class of name cflattice
#' @author Pedro Guarderas
#' @importFrom methods .valueClassTest new
#' @export
setClass( 'cflattice', slots = representation( lattice = 'list', order = 'numeric' ) )

# Extraction method --------------------------------------------------------------------------------
#' @title Extraction method
#' @description Extract method for cflattice
#' @param x object of class cflattice
#' @param i index
#' @return Numeric vector with the respective values
#' @author Pedro Guarderas
#' @export
setMethod( 
  "[", 
  "cflattice",
  function( x, i ) {
    return( x@lattice[[ i ]] )
  }
)

# Length method ------------------------------------------------------------------------------------
#' @title Length method
#' @description Extract the length of the cflattice
#' @param x object of class cflattice
#' @return Numeric value
#' @author Pedro Guarderas
#' @export
setMethod( 
  "length", 
  "cflattice",
  function( x ) {
    return( length( x@lattice ) )
  }
)

# Generic order method -----------------------------------------------------------------------------
#' @title Generic order method
#' @description Extract the order of cflattice
#' @param x object of class cflattice
#' @return Numeric value
#' @author Pedro Guarderas
#' @export
setGeneric( 
  "order", 
  valueClass = "numeric", 
  function( x ) {
    standardGeneric( "order" )
  }
)

# Order method -------------------------------------------------------------------------------------
#' @title Order method
#' @description Extract the order of cflattice
#' @param x object of class cflattice
#' @return Numeric value
#' @author Pedro Guarderas
#' @export
setMethod( 
  "order", 
  "cflattice",
  function( x ) {
    return( x@order )
  }
)
