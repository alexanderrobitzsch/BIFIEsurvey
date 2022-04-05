## File Name: svrepdesign2datalist.R
## File Version: 0.081

svrepdesign2datalist <- function(svrepdesign, varnames=NULL)
{
    if (inherits(svrepdesign,"svyimputationList")){
        datalist <- list()
        designs <- svrepdesign$designs
        Nimp <- length(designs)
        if (is.null(varnames)){
            varnames <- setdiff( colnames(designs[[1]]$variables), "one")
        }
        for (ii in 1:Nimp){
            data_ii <- designs[[ii]]$variables
            datalist[[ii]] <- data_ii[,varnames, drop=FALSE]
        }
    }
    if (inherits(svrepdesign,"svyrep.design")){
        if (is.null(varnames)){
            varnames <- setdiff( colnames(svrepdesign$variables), "one")
        }
        datalist <- list()
        data1 <- svrepdesign$variables
        datalist[[1]] <- data1[,varnames, drop=FALSE]
    }
    return(datalist)
}
