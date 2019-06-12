## File Name: BIFIE.lavaan.survey.R
## File Version: 0.641


BIFIE.lavaan.survey <- function(lavmodel, svyrepdes, lavaan_fun="sem",
    lavaan_survey_default=FALSE, fit.measures=NULL, ...)
{
    CALL <- match.call()
    s1 <- Sys.time()
    NMI <- FALSE

    #* define fit statistics
    fit.measures <- BIFIE_lavaan_survey_define_fit_measures(fit.measures=fit.measures)

    #* handle design
    is_survey_design <- FALSE
    NMI <- FALSE
    Nimp_NMI <- NULL
    variables <- NULL
    if ( class(svyrepdes)=="svyrep.design" ){
        svyrepdes0 <- svyrepdes
        data0 <- as.data.frame(svyrepdes$variables)
        Nimp <- 1
        fayfac <- svyrepdes0$scale
        lavaan_survey_default <- TRUE
        RR <- ncol(svyrepdes0$repweights)
        is_survey_design <- TRUE
    }
    if ( class(svyrepdes)=="svyimputationList"){
        svyrepdes0 <- svyrepdes$designs[[1]]
        data0 <- as.data.frame(svyrepdes0$variables)
        Nimp <- length(svyrepdes$designs)
        fayfac <- svyrepdes0$scale
        RR <- ncol(svyrepdes0$repweights)
        is_survey_design <- TRUE
    }
    if ( class(svyrepdes)=="BIFIEdata"){
        data0 <- svyrepdes$dat1
        Nimp <- svyrepdes$Nimp
        fayfac <- svyrepdes$fayfac
        NMI <- svyrepdes$NMI
        if (NMI){
            lavaan_survey_default <- FALSE
        }
        Nimp_NMI <- svyrepdes$Nimp_NMI
        svyrepdes$NMI <- FALSE
        RR <- svyrepdes$RR
        if (lavaan_survey_default){
            svyrepdes <- BIFIEdata2svrepdesign(bifieobj=svyrepdes)
        } else {
            variables <- BIFIE_lavaan_survey_define_variables(lavmodel=lavmodel,
                                svyrepdes=svyrepdes)
            datalist <- BIFIE.BIFIEdata2datalist(bifieobj=svyrepdes, varnames=variables)
        }
    }
    N <- nrow(data0)

    #- fit initial lavaan model
    lav_fun <- BIFIE_lavaan_survey_define_lavaan_function(lavaan_fun=lavaan_fun)
    lavmodel__ <- lavmodel
    args <- list(x="lavmodel__", value=lavmodel, pos=1)
    res <- do.call(what="assign", args=args)
    lavfit <- lav_fun(lavmodel__, data=data0, ...)
    class_lav <- class(lavfit)
    lavfit_coef <- BIFIE_lavaan_coef(object=lavfit)
    npar <- length(lavfit_coef)

    #* wrapper to lavaan.survey
    if (lavaan_survey_default){
        res <- BIFIE_lavaan_survey_lavaan_survey(lavaan.fit=lavfit, survey.design=svyrepdes)
        fitstat <- BIFIE_lavaan_fitMeasures(object=res, fit.measures=fit.measures)
        results <- BIFIE_lavaan_coef(object=res)
        variances <- BIFIE_lavaan_vcov(object=res)
    } else {
        results <- list()
        variances <- list()
        fitstat <- list()
        partable <- list()
        svyrepdes0 <- NULL
        for (ii in 1:Nimp){  #-- loop over imputations
            svyrepdes0 <- BIFIE_lavaan_survey_extract_dataset(svyrepdes=svyrepdes,
                                    ii=ii, variables=variables, svyrepdes0=svyrepdes0,
                                    datalist=datalist)
            res <- BIFIE_lavaan_survey_lavaan_survey(lavaan.fit=lavfit,
                            survey.design=svyrepdes0)
            results[[ii]] <- BIFIE_lavaan_coef(object=res)
            variances[[ii]] <- BIFIE_lavaan_vcov(object=res)
            fitstat[[ii]] <- BIFIE_lavaan_fitMeasures(object=res, fit.measures=fit.measures)
            partable[[ii]] <- res@ParTable
        }

        results <- bifie_extend_list_length2(x=results)
        variances <- bifie_extend_list_length2(x=variances)

        # combine fit statistics
        fitstat <- BIFIE_lavaan_survey_combine_fit_measures(fitstat=fitstat, Nimp=Nimp)

        if (! NMI){
            # inference parameters for multiply imputed datasets
            inf_res <- BIFIE_mitools_MIcombine(results=results, variances=variances)
        } else {
            # nested multiply imputed datasets
            inf_res <- BIFIE_lavaan_survey_NMIcombine(results=results,
                            variances=variances, Nimp_NMI=Nimp_NMI)
        }

        #--- include merged parameters
        res@Fit@x <- as.vector(inf_res$coefficients)
        vcov1 <- res@vcov
        vcov1$vcov <- as.matrix(inf_res$variance)
        res@vcov <- vcov1
        # combine results for lavaan parameter table
        partable <- BIFIE_lavaan_survey_combine_partable(partable=partable,
                            Nimp=Nimp, inf_res=inf_res)
        res@ParTable <- partable
    }
    #-- output
    s2 <- Sys.time()
    time <- c(s1, s2)
    res1 <- list(lavfit=res, fitstat=fitstat, CALL=CALL, time=time,
                    NMI=NMI, fayfac=fayfac, N=N, Nimp=Nimp, Nimp_NMI=Nimp_NMI, RR=RR,
                    results=results, variances=variances, partable=partable )
    class(res1) <- "BIFIE.lavaan.survey"
    return(res1)
}
