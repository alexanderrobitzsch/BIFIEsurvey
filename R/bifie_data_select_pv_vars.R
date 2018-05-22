## File Name: bifie_data_select_pv_vars.R
## File Version: 0.06

bifie_data_select_pv_vars <- function(pvpre, cn_data)
{
    # select variables with plausible values
    nc1 <- nchar( pvpre[1] )
    pv_vars <- which( substring( cn_data, 1, nc1 )==pvpre[1] )
    #-- deselect all duplicated variables
    pv_elim <- NULL
    LP <- length(pvpre)
    for (pp in 2:LP){
        pvpre_pp <- pvpre[pp]
        pv1 <- which( substring( cn_data, 1, nchar(pvpre_pp) )==pvpre_pp )
        pv_elim <- c( pv_elim, pv1 )
    }
    pv_vars <- setdiff(pv_vars, pv_elim)
    pv_vars <- gsub( pvpre[1], "", cn_data[ pv_vars ] )
    #--- output
    return(pv_vars)
}

