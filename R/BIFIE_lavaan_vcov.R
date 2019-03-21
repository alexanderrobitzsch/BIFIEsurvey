## File Name: BIFIE_lavaan_vcov.R
## File Version: 0.04

BIFIE_lavaan_vcov <- function(object)
{
    res <- BIFIE_lavaan_lavInspect(object=object, what="vcov")
    res <- as.matrix(res)
    return(res)
}
