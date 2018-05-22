## File Name: BIFIEdata2svrepdesign.R
## File Version: 0.09

BIFIEdata2svrepdesign <- function(bifieobj)
{
    CALL <- match.call()
    Nimp <- bifieobj$Nimp
    weights <- bifieobj$wgt
    repweights <- bifieobj$wgtrep
    RR <- bifieobj$RR
    scale <- bifieobj$fayfac
    rscales <- rep(1,RR)
    if (bifieobj$NMI){
        h1 <- paste0( "Nested multiply imputed datasets cannot be converted \n",
                    " into objects for the survey package.\n")
        stop(h1)
    }

    #**** create datasets
    if (Nimp==1){
        data <- as.data.frame(bifieobj$dat1)
    }
    if (Nimp>1){
        data <- BIFIE.BIFIEdata2datalist(bifieobj=bifieobj)
        data <- mitools::imputationList(data)
    }
    #*** adjust scale factor in case of finite sampling correction
    if ( length(scale) > 1){
        rscales <- scale
        scale <- 1
    }
    #*** create svrepdesign object
    svydes <- survey::svrepdesign(data=data, weights=weights, repweights=repweights,
                    type="other", scale=scale, rscales=rscales, mse=TRUE )
    svydes$call <- CALL
    return(svydes)
}
