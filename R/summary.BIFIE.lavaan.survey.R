## File Name: summary.BIFIE.lavaan.survey.R
## File Version: 0.13


summary.BIFIE.lavaan.survey <- function(object, ... )
{
    BIFIE.summary(object)
    #- lavaan summary output
    print(BIFIE_lavaan_summary(object$lavfit, ...))
    #- fit statistics
    cat("\n\nModel Fit Statistics\n")
    print(object$fit)
}
