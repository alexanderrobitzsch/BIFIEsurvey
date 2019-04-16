## File Name: BIFIE.BIFIEcdata2BIFIEdata.R
## File Version: 0.15


#--- conversion of BIFIEcdata to BIFIEdata object
BIFIE.BIFIEcdata2BIFIEdata <- function( bifieobj, varnames=NULL, impdata.index=NULL )
{
    if ( ! bifieobj$cdata ){
        stop( "You may want to use 'BIFIE.BIFIEdata2BIFIEcdata'\n")
    }

    #*** select some imputed datasets or some variables
    bifieobj <- BIFIE.cdata.select( bifieobj=bifieobj, varnames=varnames,
                    impdata.index=impdata.index )
    #*** conversion to BIFIEdata object
    bifieobj$datalistM <- bifiesurvey_rcpp_bifiecdata2bifiedata(
                                datalistM_ind=as.matrix(bifieobj$datalistM_ind),
                                datalistM_imputed=as.matrix(bifieobj$datalistM_imputed),
                                Nimp=bifieobj$Nimp, dat1=as.matrix(bifieobj$dat1),
                                datalistM_impindex=as.matrix(bifieobj$datalistM_impindex) )$datalistM
    bifieobj$cdata <- FALSE
    bifieobj$datalistM_imputed <- NULL
    bifieobj$datalistM_impindex <- NULL
    bifieobj$datalistM_ind <- NULL
    bifieobj$wgtrep <- as.matrix(bifieobj$wgtrep)
    return(bifieobj)
}
