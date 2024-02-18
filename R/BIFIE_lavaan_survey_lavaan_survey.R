## File Name: BIFIE_lavaan_survey_lavaan_survey.R
## File Version: 0.078

BIFIE_lavaan_survey_lavaan_survey <- function(lavaan.fit, survey.design, ...)
{
    args <- c(as.list(environment()), list(...))
    do.call(what=requireNamespace, args=list(package="lavaan.survey"))
    res <- do.call(what="lavaan.survey", args=args)
    return(res)
}
