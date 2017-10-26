## File Name: BIFIEdata.select.R
## File Version: 1.02
## File Last Change: 2017-09-21 17:54:23


#..............................................................
# wrapper function for subfunctions BIFIE.data.select and
# BIFIE.cdata.select
BIFIEdata.select <- function( bifieobj , varnames = NULL , impdata.index = NULL ){
	cdata <- bifieobj$cdata
	if ( cdata ){
			bifieobj <- BIFIE.cdata.select( bifieobj=bifieobj , 
				varnames=varnames  , impdata.index=impdata.index )
				}
	if ( ! cdata ){
			bifieobj <- BIFIE.data.select( bifieobj=bifieobj , 
				varnames=varnames  , impdata.index=impdata.index )
				}
	return(bifieobj)
		}
#..............................................................


