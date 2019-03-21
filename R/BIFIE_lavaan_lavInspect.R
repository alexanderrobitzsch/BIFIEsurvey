## File Name: BIFIE_lavaan_lavInspect.R
## File Version: 0.02

BIFIE_lavaan_lavInspect <- function(object, what, ...)
{
    requireNamespace("lavaan")
    res <- lavaan::lavInspect(object=object, what=what, ...)
    return(res)
}
