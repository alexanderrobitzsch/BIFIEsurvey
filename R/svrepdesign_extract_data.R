## File Name: svrepdesign_extract_data.R
## File Version: 0.10

svrepdesign_extract_data <- function(svrepdesign, varnames=NULL)
{
    wgtrep <- svrepdesign$repweights
    fayfac <- svrepdesign$scale
    wgt <- svrepdesign$pweights
    data <- svrepdesign$variables
    N <- nrow(data)
    sv_varnames <- setdiff( colnames(data), "one")
    if (is.null(varnames)){
        varnames <- sv_varnames
    }
    RR <- ncol(wgtrep)
    #--- output
    res <- list(wgt=wgt, wgtrep=wgtrep, fayfac=fayfac, varnames=varnames,
                    data=data, N=N, RR=RR)
    return(res)
}
