## File Name: BIFIE.bifiedata2bifiecdata.R
## File Version: 2.19

#######################################################
# conversion of BIFIEdata to BIFIEcdata
BIFIE.BIFIEdata2BIFIEcdata <- function( bifieobj, varnames=NULL, impdata.index=NULL )
{
    if ( bifieobj$cdata ){
        stop( "You may want to use 'BIFIE.BIFIEcdata2BIFIEdata'\n")
    }
    #******** select some imputed datasets or some variables
    bifieobj <- BIFIE.data.select( bifieobj=bifieobj, varnames=varnames,
                    impdata.index=impdata.index )

    #**** data conversion
    res1 <- bifiesurvey_rcpp_bifiedata2bifiecdata( datalistM=bifieobj$datalistM, Nimp=bifieobj$Nimp )
    bifieobj$cdata <- TRUE
    bifieobj$datalistM <- NULL
    bifieobj$datalistM_ind <- res1$datalistM_ind
    colnames(bifieobj$datalistM_ind) <- bifieobj$varnames
    bifieobj$datalistM_imputed <- res1$datalistM_imputed
    bifieobj$datalistM_impindex <- res1$datalistM_impindex
    bifieobj$time <- Sys.time()
    return(bifieobj)
}
#######################################################
