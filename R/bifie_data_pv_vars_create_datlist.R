## File Name: bifie_data_pv_vars_create_datlist.R
## File Version: 0.02
## File Last Change: 2017-06-23 19:19:12

bifie_data_pv_vars_create_datlist <- function(pvpre, pv_vars, jktype, data)
{
	dfr <- NULL
	VV <- length(pv_vars)
	for (vv in 1:VV){
		vv1 <- pv_vars[vv]
		if (jktype != "RW_PISA"){ 
			ind.vv1 <- which( substring( colnames(data) , 1 , nchar( vv1 ) ) == pv_vars[vv] )
		} else {
			varsel <- paste0( pvpre , vv1	)
			ind.vv1 <- which( colnames(data) %in% varsel )									
		}									
		Nimp <- length(ind.vv1)
		dfr2 <- data.frame( "variable" = vv1 , "var_index" = vv , "data_index" = ind.vv1 ,
						 "impdata_index"=1:Nimp ) 				 
		dfr <- rbind( dfr , dfr2 )
	}
	sel_ind <- setdiff( 1:( ncol(data) ) , dfr$data_index )
	data0 <- data[ , sel_ind ]	
	V0 <- ncol(data0)
	newvars <- seq( V0+1 , V0+VV )
	datalist <- as.list( 1:Nimp )
	for (ii in 1:Nimp ){
		dat1 <- data.frame( data0 , data[ , dfr[ dfr$impdata_index == ii  , "data_index" ] ] )
		colnames(dat1)[ newvars ] <- pv_vars
		datalist[[ii]] <- dat1 
	}  # end imputations						
	#----- output
	return(datalist)	
}
