## File Name: bifiesurvey_lavaan_fitMeasures.R
## File Version: 0.02

bifiesurvey_lavaan_fitMeasures <- function(object, fit.measures)
{
    requireNamespace("lavaan")
    res <- lavaan::fitMeasures(object=object, fit.measures=fit.measures)
    return(res)
}
