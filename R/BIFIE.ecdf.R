## File Name: BIFIE.ecdf.R
## File Version: 0.36


#######################################################################
# empirical distribution function
BIFIE.ecdf <- function( BIFIEobj, vars, breaks=NULL, quanttype=1,
    group=NULL, group_values=NULL  )
{
    #****
    s1 <- Sys.time()
    cl <- match.call()
    bifieobj <- BIFIEobj
    if (bifieobj$cdata){
        varnames <- unique( c( vars, group, "one") )
        bifieobj <- BIFIE.BIFIEcdata2BIFIEdata( bifieobj, varnames=varnames )
    }
    FF <- Nimp <- bifieobj$Nimp
    N <- bifieobj$N
    dat1 <- bifieobj$dat1
    wgt <- bifieobj$wgt
    wgtrep <- bifieobj$wgtrep
    varnames <- bifieobj$varnames
    RR <- bifieobj$RR
    datalistM <- bifieobj$datalistM
    fayfac <- bifieobj$fayfac
    se <- FALSE

    if (RR==1){ RR <- 0 }
    if ( ! se ){
        wgtrep <- matrix( wgt, ncol=1 )
        RR <- 0
                }

    vars_index <- unlist( sapply( vars, FUN=function(vv){
                        which( varnames==vv ) }, simplify=FALSE) )
    # vars values
    VV <- length(vars)

    wgt_ <- matrix( wgt, ncol=1 )
    if ( is.null( group) ){ nogroup <- TRUE } else { nogroup <- FALSE }
    cat(paste0( "|", paste0( rep("*", FF), collapse=""), "|\n" ))
    if (nogroup){
        group <- "one"
        group_values <- c(1)
    }

    #@@@@***
    group_index <- match( group, varnames )
    #@@@@***

    if ( is.null(group_values ) ){
        t1 <- bifie_table( vec=datalistM[, group_index ] )
        group_values <- sort( as.numeric( paste( names(t1) ) ))
    }

    #@@@@***
    res00 <- BIFIE_create_pseudogroup( datalistM, group, group_index, group_values )
    res00$datalistM -> datalistM
    res00$group_index -> group_index
    res00$GR -> GR
    res00$group_values -> group_values
    res00$group -> group
    #@@@@***

    if ( is.null(breaks) ){
        breaks <- as.numeric(seq( 0, 1, .01 ))
    }

    maxval <- round( max( dat1[, vars_index ], na.rm=TRUE ) * 100 )

    #**************************************************************************#
    # Rcpp call

    res <- bifie_ecdf( datalist=datalistM, wgt1=wgt_, wgtrep=wgtrep, vars_index=vars_index-1,
                fayfac=fayfac, NI=Nimp, group_index1=group_index-1, group_values=group_values,
                breaks=breaks, quanttype=quanttype, maxval=maxval )

    #--- process output
    res <- bifie_ecdf_postproc_output( res=res, group_values=group_values, breaks=breaks,
                    VV=VV, res00=res00, vars=vars, group=group )
    ecdf_ <- res$ecdf_
    stat <- res$stat

    #@@@@***
    # multiple groupings
    stat <- BIFIE_table_multiple_groupings( dfr=stat, res00=res00 )
    #@@@@***

    #*************************** OUTPUT ***************************************
    s2 <- Sys.time()
    timediff <- c( s1, s2 )
    res1 <- list( ecdf=ecdf_, stat=stat, output=res, timediff=timediff,
            N=N, Nimp=Nimp, RR=RR, fayfac=fayfac,
            NMI=BIFIEobj$NMI, Nimp_NMI=BIFIEobj$Nimp_NMI,
            CALL=cl )
    class(res1) <- "BIFIE.ecdf"
    return(res1)
}
###################################################################################

####################################################################################
# summary for BIFIE.ecdf function
summary.BIFIE.ecdf <- function( object, digits=4, ... )
{
    BIFIE.summary(object)
    cat("Empirical Distribution Function \n")
    obji <- object$ecdf
    print_object_summary( obji, digits=digits )
}
