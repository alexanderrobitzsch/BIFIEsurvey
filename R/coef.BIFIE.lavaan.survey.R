## File Name: coef.BIFIE.lavaan.survey.R
## File Version: 0.05


coef.BIFIE.lavaan.survey <- function(object, ...)
{
    coef(object$lavfit, ...)
}


vcov.BIFIE.lavaan.survey <- function(object, ...)
{
    BIFIE_lavaan_vcov(object$lavfit, ...)
}
