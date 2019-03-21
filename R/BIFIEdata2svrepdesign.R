## File Name: BIFIEdata2svrepdesign.R
## File Version: 0.26

BIFIEdata2svrepdesign <- function(bifieobj, varnames=NULL,
        impdata.index=NULL )
{
    CALL <- match.call()
    Nimp <- bifieobj$Nimp
    weights <- bifieobj$wgt
    repweights <- bifieobj$wgtrep
    RR <- bifieobj$RR
    scale <- bifieobj$fayfac
    rscales <- rep(1,RR)
    if (bifieobj$NMI){
        mess <- paste0( "Nested multiply imputed datasets cannot be converted \n",
                    " into objects for the survey package.\n")
        stop(mess)
    }

    #**** create datasets
    if (Nimp==1){
        data <- as.data.frame(bifieobj$dat1)
        if (! is.null(varnames)){
            data <- data[,varnames, drop=FALSE]
        }
    }
    if (Nimp>1){
        data <- BIFIE.BIFIEdata2datalist(bifieobj=bifieobj, varnames=varnames,
                        impdata.index=impdata.index, as_data_frame=FALSE)
        Nimp <- length(data)
        if (Nimp==1){
            data <- data[[1]]
        } else {
            data <- mitools::imputationList(datasets=data)
        }
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
