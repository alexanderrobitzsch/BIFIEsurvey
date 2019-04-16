## File Name: BIFIE_lavaan_survey_combine_partable.R
## File Version: 0.07

BIFIE_lavaan_survey_combine_partable <- function(partable, Nimp, inf_res)
{
    partable0 <- partable[[1]]
    if (Nimp>1){
        for (ii in 2L:Nimp){
            partable_ii <- partable[[ii]]
            partable0$est <- partable0$est + partable_ii$est
        }
    }
    partable0$est <- partable0$est / Nimp
    partable <- partable0
    # include results of statistical inference
    ind_free <- which(partable$free>0)
    partable$est[ ind_free ] <- as.vector(inf_res$coefficients)
    partable$se[ ind_free ] <- sqrt(diag(as.matrix(inf_res$variance)))
    return(partable)
}
