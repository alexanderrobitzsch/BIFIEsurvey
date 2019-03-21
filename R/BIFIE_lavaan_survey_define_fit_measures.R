## File Name: BIFIE_lavaan_survey_define_fit_measures.R
## File Version: 0.05

BIFIE_lavaan_survey_define_fit_measures <- function(fit.measures)
{
    if ( is.null(fit.measures) ){
        fit.measures <- c("npar","chisq", "df",
                "pvalue", "baseline.chisq", "baseline.df", "baseline.pvalue",
                "cfi","tli", "nnfi","rfi","nfi","pnfi","ifi","rni",
                "rmsea","rmr","rmr_nomean","srmr","srmr_bentler",
                "srmr_bentler_nomean","crmr","crmr_nomean","srmr_mplus",
                "srmr_mplus_nomean","gfi","agfi","pgfi","mfi","ecvi")
    }
    return(fit.measures)
}
