## File Name: BIFIE.cdata.select.R
## File Version: 1.04
## File Last Change: 2017-01-16 14:16:42


#######################################################################
# selection variables or datasets in BIFIEcdata objects
BIFIE.cdata.select <- function( bifieobj , varnames = NULL , impdata.index = NULL ){

	if ( ! bifieobj$cdata ){
		stop("Use 'BIFIE.data.select' or the general function 'BIFIEdata.select'")
						}	
					
	# retain variable "one"
	varnames0 <- bifieobj$varnames
	if ( ! is.null(varnames) ){
			varnames <- union( varnames , intersect( "one" , varnames0) )
							}						
												
    #******* do some variable checking
    if ( ! is.null(varnames) ){
#		h1 <- setdiff( varnames , colnames(bifieobj$dat1) )
		h1 <- setdiff( varnames , bifieobj$varnames )
		
		if ( length(h1) > 0 ){ 
			stop( paste0( "Following variables not in BIFIEdata object:\n  " ,
						paste0( h1 , collapse=" " ) ) )
							}
						}
					
    #******** select some imputed datasets
    if ( ! is.null(impdata.index ) ){
        # i1 <- impdata.index - 1 
		i1 <- impdata.index
        bifieobj$datalistM_imputed <- bifieobj$datalistM_imputed[ , i1 , drop=FALSE]
#        h1 <- bifieobj$datalistM_imputed[,"_imp"]               
#        bifieobj$datalistM_imputed[,"_imp"] <- match( h1 , i1 ) - 1
        bifieobj$Nimp <- length(i1)             
                }
								
    #********* select some variables
    if ( ! is.null( varnames) ){	
		
        dfr1 <- data.frame( "varnames" = bifieobj$varnames , 
					"index" = seq(1,length(bifieobj$varnames) ) )
        dfr1$selectvars <- 1 * ( dfr1$varnames %in% varnames )		
        dfr1 <- dfr1[ dfr1$selectvars == 1 , ]		
        bifieobj$datalistM_ind <- bifieobj$datalistM_ind[ , dfr1$index ]	
		i1 <- bifieobj$datalistM_impindex[,2] %in% ( dfr1$index - 1 )
		bifieobj$datalistM_imputed <- bifieobj$datalistM_imputed[ i1 , , drop=FALSE]
        bifieobj$datalistM_impindex <- bifieobj$datalistM_impindex[ i1 , , drop=FALSE]		
		
		impindex2 <- match( bifieobj$datalistM_impindex[,2] , dfr1$index - 1 ) - 1
		bifieobj$datalistM_impindex[,2] <- impindex2
				
#        bifieobj$datalistM_imputed[,"variable"] <- 
#                    match( bifieobj$datalistM_imputed[,"variable"] + 1 , dfr1$index ) - 1           
        bifieobj$dat1 <- bifieobj$dat1[ , dfr1$index , drop=FALSE]  
        bifieobj$varnames <- bifieobj$varnames[ dfr1$index ]
		
		# process variable list	
		bifieobj$variables <- bifieobj$variables[  dfr1$index , , drop=FALSE]			
    }
					
	bifieobj$Nvars <- ncol(bifieobj$dat1)					
    return(bifieobj)
        }
############################################################################
