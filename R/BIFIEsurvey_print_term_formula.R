## File Name: BIFIEsurvey_print_term_formula.R
## File Version: 0.02

BIFIEsurvey_print_term_formula <- function(formula)
{
    res <- paste0( attr( stats::terms(formula), "term.labels" ), collapse=" " )
    return(res)
}
