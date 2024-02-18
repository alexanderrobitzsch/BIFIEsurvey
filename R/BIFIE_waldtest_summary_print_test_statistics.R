## File Name: BIFIE_waldtest_summary_print_test_statistics.R
## File Version: 0.051

BIFIE_waldtest_summary_print_test_statistics <- function(object, digits,
    value_name="stat.D")
{
    if ( ! object$NMI ){ cat("D1 and D2 Statistic for Wald Test \n\n") }
    if ( object$NMI ){ cat("D1 Statistic for Wald Test \n\n") }
    obji <- object[[ value_name ]]
    print_object_summary( obji, digits=digits )
}
