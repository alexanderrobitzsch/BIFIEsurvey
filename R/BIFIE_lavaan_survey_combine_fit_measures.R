## File Name: bifie_lavaan_survey_combine_fit_measures.R
## File Version: 0.11

bifie_lavaan_survey_combine_fit_measures <- function(fitstat, Nimp)
{
    fitstat1 <- 0
    for (ii in 1:Nimp){
        fitstat1 <- fitstat1 + fitstat[[ii]]
    }
    fitstat <- fitstat1 / Nimp
    return(fitstat)
}
