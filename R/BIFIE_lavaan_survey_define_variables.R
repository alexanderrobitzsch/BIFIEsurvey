## File Name: BIFIE_lavaan_survey_define_variables.R
## File Version: 0.04

BIFIE_lavaan_survey_define_variables <- function(lavmodel, svyrepdes)
{
    requireNamespace("lavaan")
    res <- lavaan::lavaanify(model=lavmodel)
    variables0 <- paste(svyrepdes$variables$variable)
    variables <- intersect( variables0, unique( c( res$lhs, res$rhs) ) )
    return(variables)
}
