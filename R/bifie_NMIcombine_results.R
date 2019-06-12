## File Name: BIFIE_NMIcombine_results.R
## File Version: 0.09

BIFIE_NMIcombine_results <- function(results, Nimp_NMI, package="stats")
{
    if (package=="stats"){
        fun_coef <- coef
        fun_vcov <- vcov
    }
    if (package=="lavaan"){
        fun_coef <- BIFIE_lavaan_coef
        fun_vcov <- BIFIE_lavaan_vcov
    }
    #- estimates
    qhat <- BIFIE_NMIcombine_results_extract_parameters(results=results,
                    fun=fun_coef, Nimp_NMI=Nimp_NMI)

    #- variance matrices
    u <- BIFIE_NMIcombine_results_extract_parameters(results=results,
                    fun=fun_vcov, Nimp_NMI=Nimp_NMI)

    #- inference
    stat <- miceadds::NMIcombine(qhat=qhat, u=u)
    return(stat)
}
