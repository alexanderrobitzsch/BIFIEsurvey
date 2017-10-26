## File Name: bifie_table.R
## File Version: 1.07
## File Last Change: 2017-10-26 13:47:56

###########################################
# Rcpp version of R's table function
bifie_table <- function( vec  , sort.names=FALSE )
{
    datavec <- matrix( vec , ncol=1 )
    # res <- bifie_fasttable( datavec )	
	if ( storage.mode(vec) == "character" ){
				characters <- TRUE	
					} else {
				characters <- FALSE
						}
	if ( ! characters ){
		res <- bifie_fasttable( datavec )
		res1 <- res$tableM[ 1:res$N_unique , , drop=FALSE]
		tvec <- res1[,2]
		names(tvec) <- res1[,1]
				}
	if ( characters ){ 			
		t1 <- bifie_table1_character( vec )
		res <- t1$tableM
		names(res) <- t1$table_names
		if ( sort.names ){ 
			tvec <- res[ sort( names(res) ) ]
					} else { tvec <- res }
					}								
    return(tvec)
        }
###########################################		

bifietable <- bifie_table
fasttable <- bifie_table 
