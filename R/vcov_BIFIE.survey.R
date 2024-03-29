## File Name: vcov_BIFIE.survey.R
## File Version: 0.297


#--- vcov_BIFIEsurvey
vcov_BIFIEsurvey <- function( object, type=NULL, eps=1E-10, avoid.singul=FALSE )
{
    # extract replicated parameters
    parsres <- extract.replicated.pars( BIFIE.method=object, type=type )
    res1 <- object
    #*****
    parsM <- parsres$parsM
    parsrepM <- parsres$parsrepM
    parnames <- parsres$parnames
    RR <- object$RR
    if ( inherits(object,"BIFIE.correl") & is.null(type) ){
        avoid.singul <- TRUE
    }
    if ( inherits(object, c("BIFIE.freq", "BIFIE.crosstab") )  ){
        avoid.singul <- TRUE
    }
    if ( avoid.singul ){
        parsrepM <- parsrepM * ( 1 + stats::runif(prod(dim(parsrepM)), 0, eps ) )
    }
    fayfac <- res1$fayfac
    NP <- nrow(parsM)
    Cdes <- matrix( 1, ncol=NP, nrow=1 )
    Ccols <- which( colSums( abs( Cdes) ) > 0 )
    rdes <- c(0)
    # compute covariance matrices
    res0 <- bifie_comp_vcov(  parsM=parsM, parsrepM=parsrepM,
                        Cdes, rdes, Ccols - 1, fayfac=fayfac )
    var_w <- res0$var_w
    var_b <- res0$var_b
    Nimp <- res1$Nimp

    # total variance
    var_tot <- var_w  + ( 1 + 1/Nimp ) * var_b
    rownames(var_tot) <- colnames(var_tot) <- parnames
    if (object$NMI){
        var_tot <- BIFIE_NMI_inference_parameters( parsM, parsrepM, fayfac,
                RR, Nimp, object$Nimp_NMI, comp_cov=TRUE )$Tm
    }
    return(var_tot)
}

vcov.BIFIE.correl <- function( object, type=NULL, ... )
{
    pars <- vcov_BIFIEsurvey( object=object, type=type, ... )
    return(pars)
}
# further BIFIE functions
vcov.BIFIE.by <- function( object, ... )
{
    pars <- vcov_BIFIEsurvey( object=object, type=NULL, ...)
    return(pars)
}
vcov.BIFIE.derivedParameters <- function( object, ... )
{
    pars <- vcov_BIFIEsurvey( object=object, type=NULL, ... )
    return(pars)
}
vcov.BIFIE.crosstab <- function( object, ... )
{
    pars <- vcov_BIFIEsurvey( object=object, type=NULL, ... )
    return(pars)
}
vcov.BIFIE.freq <- function( object, ... )
{
    pars <- vcov_BIFIEsurvey( object=object, type=NULL, ...)
    return(pars)
}
vcov.BIFIE.linreg <- function( object, ... )
{
    pars <- vcov_BIFIEsurvey( object=object, type=NULL, ... )
    return(pars)
}
vcov.BIFIE.logistreg <- function( object, ... )
{
    pars <- vcov_BIFIEsurvey( object=object, type=NULL, ...)
    return(pars)
}
vcov.BIFIE.univar <- function( object, ... )
{
    pars <- vcov_BIFIEsurvey( object=object, type=NULL, ...)
    return(pars)
}
vcov.BIFIE.twolevelreg <- function( object, ... )
{
    if (object$se){
        pars <- vcov_BIFIEsurvey( object=object, type=NULL, ... )
    } else {
        pars <- vcov( object$micombs )
    }
    return(pars)
}
vcov.BIFIE.pathmodel <- function( object, ... )
{
    pars <- vcov_BIFIEsurvey( object=object, type=NULL, ...)
    return(pars)
}

