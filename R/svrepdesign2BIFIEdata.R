## File Name: svrepdesign2BIFIEdata.R
## File Version: 0.261


svrepdesign2BIFIEdata <- function(svrepdesign, varnames=NULL, cdata=FALSE)
{
    ## class svyrep.design
    if (inherits(svrepdesign,"svyrep.design")){
        res <- svrepdesign_extract_data(svrepdesign=svrepdesign, varnames=varnames)
        wgt <- res$wgt
        wgtrep <- res$wgtrep
        fayfac <- res$fayfac
        varnames <- res$varnames
        data <- svrepdesign$variables
        data$one <- NULL
        datalist <- data[, varnames, drop=FALSE]
    }
    ## class svyimputationList
    if (inherits(svrepdesign,"svyimputationList")){
        designs <- svrepdesign$designs
        Nimp <- length(designs)
        svrepdesign0 <- designs[[1]]
        res <- svrepdesign_extract_data(svrepdesign=svrepdesign0, varnames=varnames)
        wgt <- res$wgt
        wgtrep <- res$wgtrep
        fayfac <- res$fayfac
        varnames <- res$varnames
        datalist <- svrepdesign2datalist(svrepdesign=svrepdesign, varnames=varnames)
    }
    #- convert to BIFIEdata object
    res <- BIFIE.data(data.list=datalist, wgt=wgt, wgtrep=wgtrep, fayfac=fayfac,
                    cdata=cdata, NMI=FALSE)
    return(res)
}
