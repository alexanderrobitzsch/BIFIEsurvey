## File Name: BIFIE_lavaan_fitMeasures.R
## File Version: 0.03

BIFIE_lavaan_fitMeasures <- function(object, fit.measures)
{
    requireNamespace("lavaan")
    res <- lavaan::fitMeasures(object=object, fit.measures=fit.measures)
    return(res)
}
