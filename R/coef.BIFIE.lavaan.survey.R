## File Name: coef.BIFIE.lavaan.survey.R
## File Version: 0.03


coef.BIFIE.lavaan.survey <- function(object, ...)
{
    coef(object$lavfit, ...)
}


vcov.BIFIE.lavaan.survey <- function(object, ...)
{
    vcov(object$lavfit, ...)
}
