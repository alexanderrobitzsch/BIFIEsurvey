## File Name: bifie_ecdf_postproc_output.R
## File Version: 0.11


bifie_ecdf_postproc_output <- function(res, group_values, breaks, VV, res00,
    vars, group)
{
    GG <- length(group_values)
    BB <- length(breaks)
    #--- output data frame containing distribution functions as columns
    dfr <- as.data.frame( matrix( NA, nrow=BB, ncol=VV*GG+1 ) )
    colnames(dfr)[1] <- "yval"
    dfr[,1] <- breaks
    dfr[,-1] <- matrix( res$ecdf, nrow=BB, ncol=VV*GG )
    colnames(dfr)[-1] <- paste0( rep( vars, each=GG ),
                        "_", group, rep( group_values, VV ) )
    GR <- res00$GR
    if (GR>1){
        group_orig <- res00$group_orig
        group_values_recode <- res00$group_values_recode
        rr <- 1
        p1 <- paste0( group_orig[rr], group_values_recode[,rr] )
        for (rr in 2:GR){
            p1 <- paste0( p1, "_", group_orig[rr], group_values_recode[,rr] )
        }
        cn <- paste0( rep( vars, each=GG ), "_", rep( p1, VV ) )
        colnames(dfr)[-1] <- cn
    }
    ecdf_ <- dfr
    #--- data frame with statistics
    stat <- NULL
    ii <- 1
    for (vv in 1:VV){
        for (gg in 1:GG){
            dfr1 <- data.frame( "var"=vars[vv], "groupvar"=group,
                    "groupval"=group_values[gg], "yval"=breaks, "quant"=dfr[,ii+1] )
            ii <- ii + 1
            stat <- rbind( stat, dfr1 )
        }
    }
    #--- output
    res <- list( ecdf_=ecdf_, stat=stat)
    return(res)
}

