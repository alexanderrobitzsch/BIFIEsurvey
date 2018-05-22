## File Name: load.BIFIEdata.files.R
## File Version: 3.17


#****** loading files for conversion into BIFIE data objects
load.BIFIEdata.files <- function( files.imp, wgt, file.wgtrep, file.ind=NULL, type="Rdata",
        varnames=NULL, cdata=TRUE, dir=getwd(), ... )
{
    Nimp <- length(files.imp)
    # handle cases in which no weights are defined?

    #**** no indicator dataset
    if ( is.null(file.ind) ){
        #**********
        # read imputed datasets
        datalist <- list(1:Nimp)
        for (ii in 1:Nimp){
            # cat(paste0( "- Read ", files.imp[[ii]], "\n") )
            # utils::flush.console()
            load_BIFIEdata_files_cat_print_file_name( file=files.imp[[ii]] )
            dat1 <- miceadds::load.data( file=files.imp[[ii]], path=dir, type=type, ... )
            if (ii==1){  wgt <- dat1[, wgt ] }
            dat1 <- as.data.frame(dat1)
            dat1 <- load_BIFIEdata_files_select_variables( dat=dat1, varnames=varnames )
            dat1$one <- NULL
            datalist[[ii]] <- dat1
        }
        #**************
        # read replicate weights
        load_BIFIEdata_files_cat_print_file_name( file=file.wgtrep )
        wgtrep <- miceadds::load.data( file=file.wgtrep, type=type, path=dir, ...)
        #****************************
        # create BIFIEdata object
        bifieobj <- BIFIE.data( data.list=datalist, wgt=wgt, wgtrep=wgtrep, cdata=cdata )
    }

    #**** with indicator dataset
    if ( ! is.null( file.ind ) ){
        #--- read indicator dataset
        load_BIFIEdata_files_cat_print_file_name( file=file.ind )
        dat_ind <- miceadds::load.data( file=file.ind, type=type, path=dir, ...)
        if ( is.null(varnames) ){
            varnames <- setdiff( colnames(dat_ind ), "one" )
        }
        dat_ind <- load_BIFIEdata_files_select_variables( dat=dat_ind, varnames=varnames )
        dat_ind <- as.matrix( dat_ind )
        # add column 1 for "one"
        dat_ind <- cbind( dat_ind, 1 )
        colnames(dat_ind) <- c( varnames, "one")

        #************************
        # Read first imputed dataset
        ii <- 1
        load_BIFIEdata_files_cat_print_file_name( file=files.imp[[ii]] )
        dat1 <- miceadds::load.data( file=files.imp[[ii]], path=dir, type=type, ... )
        dat1 <- load_BIFIEdata_files_select_variables( dat=dat1, varnames=varnames )
        datalist <- list( dat1 )
        #**************
        # read replicate weights
        load_BIFIEdata_files_cat_print_file_name( file=file.wgtrep )
        wgtrep <- miceadds::load.data( file=file.wgtrep, type=type, path=dir, ...)
        #***************
        # create initial BIFIEdata object
        bifieobj <- BIFIE.data( data.list=datalist, wgt=wgt, wgtrep=wgtrep, cdata=cdata )
        bifieobj$datalistM_ind <- dat_ind
        bifieobj$Nimp <- Nimp
        Nmiss <- sum( 1 - dat_ind )
        datalistM_imputed <- matrix( NA, nrow=Nmiss, Nimp)
        res1 <- bifiesurvey_rcpp_bifiedata_stepwise( as.matrix(dat1), dat_ind, Nmiss )$datalistM_imputed
        datalistM_imputed[,1] <- res1[,4]
        datalistM_impindex <- res1[,2:3]

        #----- read other imputed datasets
        if (Nimp>1){
            for (ii in 2:Nimp){
                load_BIFIEdata_files_cat_print_file_name( file=files.imp[[ii]] )
                dat1 <- miceadds::load.data( file=files.imp[[ii]], path=dir, type=type, ... )
                dat1 <- load_BIFIEdata_files_select_variables( dat=dat1, varnames=varnames )
                res1 <- bifiesurvey_rcpp_bifiedata_stepwise(  as.matrix(dat1), dat_ind, Nmiss )$datalistM_imputed
                datalistM_imputed[,ii] <- res1[,4]
                datalistM_impindex <- rbind( datalistM_impindex, res1[,2:3] )
            }
        }
        bifieobj$dat1 <- cbind( as.data.frame(dat1), "one"=1 )
        bifieobj$datalistM_imputed <- datalistM_imputed
        datalistM_imputed <- NULL
        bifieobj$datalistM_impindex <- datalistM_impindex
    }  # end indicator data
    #--- return BIFIE object
    return(bifieobj)
}
