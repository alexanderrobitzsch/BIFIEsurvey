## File Name: BIFIE_lavaan_coef.R
## File Version: 0.12

BIFIE_lavaan_coef <- function(object, ...)
{
    requireNamespace("lavaan")
    lavaan_coef <- methods::getMethod("coef", "lavaan")
    est <- lavaan_coef(object, ...)
    return(est)
}


