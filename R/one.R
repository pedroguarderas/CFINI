# 
# 
# library( data.table )
# library( lubridate )
# library( Risko )
# 
# load( paste( get.data.server(), 'TablasMortalidad/RData/Tabla_Mortalidad_IESS_2000.RData', sep = '' ) )
# 
# Data<-data.table( id = 1, birth = dmy( '03-09-1989' ) )
# t0<-dmy( '01-01-2017' )
# t<-seq( dmy( '01-01-2017' ), dmy( '01-01-2040' ), by = 'year' )
# 
# daycount<-function( start, end ) {
#   d<-interval( start, end ) / dyears(1)
#   return( d )
# }
# 
# Survival<-function( Data, LifeTable, t0, t ) {
#   D<-copy( Data )
#   D[ , x := daycount( birth, t0 ) ]
#   D<-D[ x > 0 ]
#   D<-D[ data.table( t ) ]
#   D<-merge( D, data.table( t ), by = NULL, allow.cartesian = TRUE )
# }
# 
# MONT[ , X := floor( x ) ]
# MONT[ , X1 := floor( x ) + 1 ]
# MONT[ , u := x - X ]
# MONT[ , Xt := X + floor( t + u ) ]
# MONT[ , Xt1 := X + floor( t + u ) + 1 ]
# MONT[ , r := t + u - floor( t + u ) ]
# 
# MONT<-merge( MONT, tabla.mortalidad.iess[ , list( X = x, sexo, lx ) ], 
#              by = c( 'X', 'sexo' ) )
# 
# MONT<-merge( MONT, tabla.mortalidad.iess[ , list( X1 = x, sexo, lx1 = lx  ) ], 
#              by = c( 'X1', 'sexo' ) )
# 
# MONT<-merge( MONT, tabla.mortalidad.iess[ , list( Xt = x, sexo, lxt = lx ) ], 
#              by = c( 'Xt', 'sexo' ) )
# 
# MONT<-merge( MONT, tabla.mortalidad.iess[ , list( Xt1 = x, sexo, lxt1 = lx  ) ], 
#              by = c( 'Xt1', 'sexo' ) )
# 
# MONT[ lx != 0, uPx := ( lx^(1-u) ) * ( lx1^u ) / lx ]
# MONT[ lx == 0, uPx := 0 ]
# MONT[ lx != 0, tuPx := ( lxt^(1-r) ) * ( lxt1^r ) / lx ]
# MONT[ lx == 0, tuPx := 0 ]
# MONT[ uPx != 0, tPx := tuPx / uPx ]
# MONT[ uPx == 0, tPx := 0 ]