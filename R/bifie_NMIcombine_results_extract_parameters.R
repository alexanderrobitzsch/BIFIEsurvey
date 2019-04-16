## File Name: bifie_NMIcombine_results_extract_parameters.R
## File Version: 0.02


bifie_NMIcombine_results_extract_parameters <- function(results, fun, Nimp_NMI,
    loop_within=TRUE)
{
    u <- lapply(results, FUN=function(res){ fun(res) } )
    res <- miceadds::List2nestedList(List=u, N_between=Nimp_NMI[1],
                        N_within=Nimp_NMI[2], loop_within=loop_within)
    return(res)
}
