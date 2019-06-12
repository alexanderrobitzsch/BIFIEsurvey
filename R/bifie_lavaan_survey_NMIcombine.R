## File Name: BIFIE_lavaan_survey_NMIcombine.R
## File Version: 0.03

BIFIE_lavaan_survey_NMIcombine <- function(results, variances, Nimp_NMI)
{
    qhat <- miceadds::List2nestedList(List=results, N_between=Nimp_NMI[1],
                        N_within=Nimp_NMI[2], loop_within=TRUE)
    u <- miceadds::List2nestedList(List=variances, N_between=Nimp_NMI[1],
                        N_within=Nimp_NMI[2], loop_within=TRUE)
    #- inference
    inf_res <- miceadds::NMIcombine(qhat=qhat, u=u)
    inf_res$coefficients <- inf_res$qbar
    inf_res$variance <- inf_res$Tm
    return(inf_res)
}
