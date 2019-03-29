## File Name: BIFIE_lavaan_survey_lavaan_survey.R
## File Version: 0.05

BIFIE_lavaan_survey_lavaan_survey <- function(lavaan.fit, survey.design, ...)
{
    requireNamespace("lavaan.survey")
    res <- lavaan.survey::lavaan.survey(lavaan.fit=lavaan.fit,
                survey.design=survey.design, ...)
    return(res)
}
