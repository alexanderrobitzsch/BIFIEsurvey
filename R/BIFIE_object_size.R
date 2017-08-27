## File Name: BIFIE_object_size.R
## File Version: 0.01
## File Last Change: 2017-06-23 19:46:24


##################################################
# output object sizes in MB and GB
BIFIE_object_size <- function( x1 ){
	x1 <- utils::object.size(x1)
	vals <- c( x1 , as.numeric(x1 / 1024^2) ,
					as.numeric(x1 / 1024^3)  )
	names(vals) <- c("B" , "MB" , "GB" )		
    res <- list( "value" = vals )
	res$value_string <- paste0( round( vals[2] , 3 ) , 
			" MB | " , round( vals[3] , 5 ) ,	" GB" ) 		
    return(res)
}
