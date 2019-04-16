## File Name: BIFIE.BIFIEdata2datalist.R
## File Version: 0.19


#--- converts a BIFIEdata object into a list of multiply imputed datasets
BIFIE.BIFIEdata2datalist <- function( bifieobj, varnames=NULL,
        impdata.index=NULL, as_data_frame=FALSE )
{
    Nimp <- bifieobj$Nimp
    NMI <- bifieobj$NMI
    bifieobj <- BIFIEdata.select(bifieobj=bifieobj, varnames=varnames,
                        impdata.index=impdata.index )
    if (bifieobj$cdata){
        bifieobj <- BIFIE.BIFIEcdata2BIFIEdata( bifieobj=bifieobj)
    }
    datalistM <- bifieobj$datalistM
    variables <- bifieobj$variables
    cndat1 <- colnames(bifieobj$dat1)
    N <- bifieobj$N
    Nimp <- bifieobj$Nimp
    datalist <- as.list(1:Nimp)
    for (ii in 1:Nimp){
        dat0 <- datalistM[ (ii-1)*N + 1:N, ]
        colnames(dat0) <- cndat1
        datalist[[ii]] <- as.data.frame(dat0)
    }
    if ((Nimp==1) & as_data_frame){
        datalist <- datalist[[1]]
    }
    return(datalist)
}

