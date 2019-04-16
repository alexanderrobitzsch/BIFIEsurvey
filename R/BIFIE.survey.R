## File Name: BIFIE.survey.R
## File Version: 0.224

BIFIE.survey <- function(svyrepdes, survey.function, ...)
{
    CALL <- match.call()
    s1 <- Sys.time()
    NMI <- FALSE
    svrepdes <- svyrepdes
    if ( class(svyrepdes)=="BIFIEdata"){
        data0 <- svyrepdes$dat1
        N <- nrow(data0)
        Nimp <- svyrepdes$Nimp
        fayfac <- svyrepdes$fayfac
        NMI <- svyrepdes$NMI
        Nimp_NMI <- svyrepdes$Nimp_NMI
        svyrepdes$NMI <- FALSE
        RR <- svyrepdes$RR
        wgt <- svyrepdes$wgt
        wgtrep <- svyrepdes$wgtrep
        variables <- NULL
        args <- list(...)
        for (vv in c("formula", "x")){
            if ( vv %in% names(args)){
                args_vv <- args[[vv]]
                if (class(args_vv)=="formula"){
                    variables <- all.vars(args_vv)
                }
            }
        }
        datalist <- BIFIE.BIFIEdata2datalist( bifieobj=svyrepdes, varnames=variables)
    }
    if ( class(svyrepdes)=="svyimputationList"){
        res <- svrepdesign_extract_data(svrepdesign=svrepdes$designs[[1]])
        N <- res$N
        RR <- res$RR
        fayfac <- res$fayfac
        Nimp <- length(svrepdes$designs)
    }

    #* loop over imputations
    if ( class(svyrepdes) %in% c("BIFIEdata", "svyimputationList") ){
        res <- list()
        svyrep_ii <- NULL
        for (ii in 1:Nimp){
            if ( class(svyrepdes)=="BIFIEdata"){
                svyrep_ii <- BIFIE_lavaan_survey_extract_dataset( svyrepdes=svyrepdes,
                                    ii=ii, variables=NULL, svyrepdes0=svyrep_ii, datalist=datalist)
            }
            if ( class(svyrepdes)=="svyimputationList"){
                svyrep_ii <- svrepdes$designs[[ii]]
            }
            args <- list(...)
            args$design <- svyrep_ii
            res[[ii]] <- do.call( what=survey.function, args=args)
        }
        results <- res
        results <- bifie_extend_list_length2(x=results)
    }

    if (! NMI){
        #*** statistical inference using mitools package
        stat <- BIFIE_mitools_MIcombine(results=results)
    } else {
        #*** nested multiply imputed dataset
        stat <- bifie_NMIcombine_results(results=results, Nimp_NMI=Nimp_NMI, package="stats")
    }

    #-- output
    s2 <- Sys.time()
    time <- c(s1, s2)
    res1 <- list(stat=stat, CALL=CALL, time=time, NMI=NMI, fayfac=fayfac, N=N,
                Nimp=Nimp, RR=RR, results=results, Nimp_NMI=Nimp_NMI)
    class(res1) <- "BIFIE.survey"
    return(res1)
}


#-- summary function
summary.BIFIE.survey <- function( object, digits=3, ... )
{
    BIFIE.summary(object)
    cat("Estimated Parameters\n")
    summary(object$stat, digits=digits)
}
