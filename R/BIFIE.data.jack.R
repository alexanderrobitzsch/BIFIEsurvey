## File Name: BIFIE.data.jack.R
## File Version: 1.67
###########################################################
# BIFIE.data objects for designs with jackknife zones
BIFIE.data.jack <- function( data, wgt=NULL, jktype="JK_TIMSS", pv_vars=NULL,
    jkzone=NULL, jkrep=NULL, jkfac=NULL, fayfac=NULL,
    wgtrep="W_FSTR", pvpre=paste0("PV",1:5), ngr=100,
    seed=.Random.seed,
    cdata=FALSE )
{
    cl <- match.call()

    # subroutine for preparation of nested multiple imputations
    # res0 <- BIFIE_data_nested_MI( data.list=data.list, NMI=NMI )
    # data.list <- res0$data.list
    # Nimp_NMI <- res0$Nimp_NMI

    fayfac0 <- fayfac

    if ( ( ! is.null(wgtrep) ) & ( is.null(fayfac) ) ){
        fayfac <- 1
    }

    #*** list of multiply imputed datasets
    if ( ( is.list(data) ) & ( ! is.data.frame(data) ) ){
        dataL <- data
        data <- dataL[[1]]
    } else {
        dataL <- data
    }
    data <- as.data.frame( data )
    #*********************************************************
    # using fixed jackknife zones
    if (jktype=="JK_GROUP"){
        N <- nrow(data)
        if ( is.null(wgt) ){
            data$wgt <- rep(1,N)
            wgt <- "wgt"
        }
        data$jkrep <- rep(0,N)
        jkrep <- "jkrep"
        fayfac <- ngr / ( ngr - 1 )
        jkfac <- 0
    }


    #**********************************************************
    #*** defaults for jackknife creation: random groups
    if (jktype=="JK_RANDOM"){
        N <- nrow(data)
        if ( is.null(wgt) ){
            data$wgt <- rep(1,N)
            wgt <- "wgt"
        }
        if ( ! is.null(seed) ){
            set.seed( seed )
            indzone <- sample(1:N)
        } else {
            indzone <- 1:N
        }
        jkzone <- 1:N
        N1 <- N / ngr
        jkzone <- floor( jkzone / ( N1 + 1E-5  ) ) + 1
        jkzone <- jkzone[indzone]
        jkrep <- rep(0,N)
        data$jkzone <- jkzone
        jkzone <- "jkzone"
        data$jkrep <- jkrep
        jkrep <- "jkrep"
        fayfac <- ngr / ( ngr - 1 )
        jkfac <- 0
    }

    #**********************************************************
    #**** defaults for TIMSS
    if (jktype %in% c("JK_TIMSS","JK_TIMSS2") ){
        if ( is.null(jkrep) ){
            jkrep <- "JKREP"
        }
        if ( is.null(jkzone) ){
            jkzone <- "JKZONE"
        }
        if ( is.null(wgt) ){
            wgt <- "TOTWGT"
        }
        jkfac <- 2
    }
    #***********************************************************
    #**** defaults for PISA
    if (jktype=="RW_PISA"){
        jkrep <- NULL
        jkzone <- NULL
        if ( is.null(wgt)){
            wgt <- "W_FSTUWT"
        }
        jkfac <- NULL
        cn_data <- colnames(data)
        repvars <- grep( wgtrep, cn_data )
        RR <- length(repvars)

        pv_vars <- bifie_data_select_pv_vars(pvpre, cn_data )
        datarep <- data[, repvars ]
        RR <- ncol(datarep)
        fayfac <- 1 /  RR / ( 1 - .5)^2
        data <- data[, - repvars ]
    }
    #******** generate replicate weights
    if ( jktype %in% c("JK_TIMSS", "JK_GROUP", "JK_RANDOM", "JK_TIMSS2") ) {
        # redefine jackknife zones
        jkzones1 <- unique( data[,jkzone] )
        data[,jkzone] <- match( data[,jkzone], jkzones1)
        #***********
        RR <- max( data[,jkzone] )
        prblen <- 10
        prbar <- BIFIE.progressbar( ops=RR, prblen=prblen )
        cat("+++ Generate replicate weights\n")
        cat(paste0("|", paste0(rep("*",prblen), collapse=""), "|\n|"))
        utils::flush.console()
        addname <- 10^( floor( log( RR+.5, 10 ) )  + 1 )
        data[, jkzone ] <- match( data[, jkzone ], unique( data[, jkzone] ) )
        datarep <- bifiesurvey_rcpp_jackknife_timss( wgt=data[,wgt], jkzone=data[,jkzone]-1,
                        jkrep=data[,jkrep], RR=RR, jkfac=jkfac,  prbar=prbar )
        colnames(datarep) <- paste0("w_fstr", substring( paste0(addname +1:RR),2) )
        # adjustments for JK_TIMSS2 type
        if (jktype=="JK_TIMSS2"){
            datarep0 <- bifiesurvey_rcpp_jackknife_timss( wgt=data[,wgt], jkzone=data[,jkzone]-1,
                            jkrep=1 - data[,jkrep], RR=RR, jkfac=jkfac,  prbar=prbar )
            colnames(datarep0) <- paste0("w_fstr", substring( paste0(addname +1:RR),2) )
            datarep <- cbind( datarep, datarep0 )
            ind_rep <- unlist( sapply( 1:RR, FUN=function(rr){ rr + c(0,RR) },
                                simplify=FALSE) )
            datarep <- datarep[, ind_rep ]
            addname <- 10^( floor( log( 2*RR+.5, 10 ) )  + 1 )
            colnames(datarep) <- paste0("w_fstr", substring( paste0(addname +1:(2*RR)),2) )
            RR <- 2*RR
            fayfac <- .5
        }
        cat("|\n")
    }

    #******** generate replicated datasets for datasets
    if ( is.null( pv_vars) ){
        datalist <- dataL
    }

    #--------------------------------------------------
    if ( ! is.null( pv_vars )){
        datalist <- bifie_data_pv_vars_create_datlist( pvpre=pvpre, pv_vars=pv_vars,
                            jktype=jktype, data=data )
    }  # end pv_vars
    #--------------------------------------------------

    if ( ! is.null(fayfac0) ){
        fayfac <- fayfac0
    }

    #*** create BIFIE.data object
    bifiedat <- BIFIE.data( datalist, wgt=data[, wgt ], wgtrep=datarep, fayfac=fayfac,
                            cdata=cdata, NMI=FALSE )
    bifiedat$CALL <- cl
    return(bifiedat)
}
###############################################################################
