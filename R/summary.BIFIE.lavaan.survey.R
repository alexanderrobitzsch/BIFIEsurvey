## File Name: summary.BIFIE.lavaan.survey.R
## File Version: 0.05


summary.BIFIE.lavaan.survey <- function(object, ... )
{
    BIFIE.summary(object)
    print(summary(object$lavfit, ...))
    cat("\n\nModel Fit Statistics\n")
    print(object$fit)
}
