## File Name: BIFIE_data_transform_process_formula.R
## File Version: 0.01

BIFIE_data_transform_process_formula <- function(transform.formula)
{
    t1 <- stats::terms(transform.formula)
    t2 <- attr(t1, "term.labels")
    transform.formula <- stats::as.formula( paste0( "~ 0 + ", paste0( t2, collapse=" + ") ) )
    return(transform.formula)
}
