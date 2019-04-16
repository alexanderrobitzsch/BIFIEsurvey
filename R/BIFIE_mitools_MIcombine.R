## File Name: BIFIE_mitools_MIcombine.R
## File Version: 0.02

BIFIE_mitools_MIcombine <- function(...)
{
    requireNamespace("mitools")
    res <- mitools::MIcombine(...)
    return(res)
}
