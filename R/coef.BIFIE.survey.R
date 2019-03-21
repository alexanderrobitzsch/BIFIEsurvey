## File Name: coef.BIFIE.survey.R
## File Version: 0.01


coef.BIFIE.survey <- function(object, ...)
{
    coef(object$stat, ...)
}


vcov.BIFIE.survey <- function(object, ...)
{
    vcov(object$stat, ...)
}
