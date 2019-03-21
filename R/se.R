## File Name: se.R
## File Version: 0.02

se <- function(object)
{
    res <- sqrt(diag(vcov(object)))
    return(res)
}
