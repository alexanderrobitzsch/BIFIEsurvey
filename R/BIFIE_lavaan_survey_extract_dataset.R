## File Name: BIFIE_lavaan_survey_extract_dataset.R
## File Version: 0.091

BIFIE_lavaan_survey_extract_dataset <- function(svyrepdes, ii, variables,
    svyrepdes0=NULL, datalist=NULL)
{
    if ( inherits(svyrepdes,"svyimputationList") ){
        svyrepdes0 <- svyrepdes$designs[[ii]]
    }
    if ( inherits(svyrepdes,"BIFIEdata") ){
        use_datalist <- (ii>1) & ( ! is.null(datalist) )
        if (! use_datalist){
            svyrepdes0 <- BIFIEdata2svrepdesign(bifieobj=svyrepdes, varnames=variables,
                                impdata.index=ii)
        } else {
            svyrepdes0$variables <- datalist[[ii]]
        }
    }
    return(svyrepdes0)
}
