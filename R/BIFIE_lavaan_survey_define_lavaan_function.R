## File Name: BIFIE_lavaan_survey_define_lavaan_function.R
## File Version: 0.04


BIFIE_lavaan_survey_define_lavaan_function <- function(lavaan_fun)
{
    requireNamespace("lavaan")
    if (lavaan_fun=="sem"){ lav_fun <- lavaan::sem }
    if (lavaan_fun=="cfa"){ lav_fun <- lavaan::cfa }
    if (lavaan_fun=="lavaan"){ lav_fun <- lavaan::lavaan }
    if (lavaan_fun=="growth"){ lav_fun <- lavaan::growth }
    return(lav_fun)
}
