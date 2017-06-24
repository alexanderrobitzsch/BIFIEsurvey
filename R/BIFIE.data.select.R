


#######################################################################
# selection variables or datasets in BIFIEcdata objects
BIFIE.data.select <- function( bifieobj , varnames = NULL , impdata.index = NULL ){
	if ( bifieobj$cdata ){
		stop("Use 'BIFIE.cdata.select' or the general function 'BIFIEdata.select'")
						}
	# retain variable "one"
	varnames0 <- bifieobj$varnames
	if ( ! is.null(varnames) ){
			varnames <- union( varnames , intersect( "one" , varnames0) )
							}
							
							
    #******** select some imputed datasets	
    if ( ! is.null(impdata.index ) ){
        i1 <- impdata.index - 1
		N <- bifieobj$N
		ind <- unlist( sapply( i1 , FUN = function(ii){ 
				vec <- ii*N + ( 1:N )
				return(vec)
						} , simplify=FALSE) )
		bifieobj$datalistM <- bifieobj$datalistM[ ind , , drop=FALSE]
        bifieobj$Nimp <- length(i1)             
                }
    #********* select some variables
	if ( ! is.null( varnames) ){	
		dfr1 <- data.frame( "varnames" = bifieobj$varnames , "index" = seq(1,length(bifieobj$varnames) ) )
		dfr1$selectvars <- 1 * ( dfr1$varnames %in% varnames )
		dfr1 <- dfr1[ dfr1$selectvars == 1 , ]
		bifieobj$datalistM <- bifieobj$datalistM[ , dfr1$index , drop=FALSE]
		bifieobj$dat1 <- bifieobj$dat1[ , dfr1$index , drop=FALSE]	
		bifieobj$varnames <- bifieobj$varnames[ dfr1$index ]
		# process variable list
		bifieobj$variables <- bifieobj$variables[  dfr1$index , ]						
					}
	bifieobj$Nvars <- ncol(bifieobj$dat1)
	return(bifieobj)
        }
############################################################################

