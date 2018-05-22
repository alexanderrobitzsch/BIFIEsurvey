## File Name: load_BIFIEdata_files_select_variables.R
## File Version: 0.03


load_BIFIEdata_files_select_variables <- function( dat, varnames )
{
    if ( ! is.null(varnames) ){
        dat <- dat[, varnames ]
    }
    return(dat)
}
