## File Name: BIFIE.lavaan.survey.R
## File Version: 0.46


BIFIE.lavaan.survey <- function(lavmodel, svyrepdes, lavaan_fun="sem",
    lavaan_survey_default=FALSE, ...)
{
    CALL <- match.call()
    s1 <- Sys.time()
    NMI <- FALSE
    requireNamespace("lavaan")
    requireNamespace("lavaan.survey")
    #* fit initial model with lavaan
    if ( class(svyrepdes)=="svyrep.design" ){
        svyrepdes0 <- svyrepdes
        data0 <- as.data.frame(svyrepdes$variables)
        Nimp <- 1
        lavaan_survey_default <- TRUE
    }
    if ( class(svyrepdes)=="svyimputationList"){
        svyrepdes0 <- svyrepdes$designs[[1]]
        data0 <- as.data.frame(svyrepdes0$variables)
        Nimp <- length(svyrepdes$designs)
        fayfac <- svyrepdes0$scale
    }

    # details about survey object
    N <- nrow(data0)
    fayfac <- svyrepdes0$scale
    RR <- ncol(svyrepdes0$repweights)

    #- initial lavaan model
    lav_fun <- create_function_name(pack="lavaan", fun=lavaan_fun)
    lavfit <- lav_fun(lavmodel, data=data0, ...)
    class_lav <- class(lavfit)

    #* wrapper to lavaan.survey
    if (lavaan_survey_default){
        res <- lavaan.survey::lavaan.survey(lavaan.fit=lavfit, survey.design=svyrepdes )
        fitstat <- lavaan::fitMeasures(res)
    } else {
        npar <- length(coef(lavfit))
        results <- list()
        variances <- list()
        fitstat <- 0
        for (ii in 1:Nimp){
            svyrepdes0 <- svyrepdes$designs[[ii]]
            res <- lavaan.survey::lavaan.survey(lavaan.fit=lavfit, survey.design=svyrepdes0 )
            results[[ii]] <- coef(res)
            variances[[ii]] <- vcov(res)
            fitstat <- fitstat + lavaan::fitMeasures(res)
        }
        fitstat <- fitstat / Nimp
        inf_res <- mitools::MIcombine(results=results, variances=variances)
        #--- include merged parameters
        res@Fit@x <- as.vector(inf_res$coefficients)
        res@vcov$vcov <- as.matrix(inf_res$variance)
        partable <- res@ParTable
        ind_free <- which(partable$free>0)
        partable$est[ ind_free ] <- as.vector(inf_res$coefficients)
        partable$se[ ind_free ] <- sqrt(diag(as.matrix(inf_res$variance)))
        res@ParTable <- partable
    }
    #-- output
    s2 <- Sys.time()
    time <- c(s1, s2)
    res1 <- list(lavfit=res, fitstat=fitstat, CALL=CALL, time=time,
                    NMI=NMI, fayfac=fayfac, N=N, Nimp=Nimp, RR=RR)
    class(res1) <- "BIFIE.lavaan.survey"
    return(res1)
}
