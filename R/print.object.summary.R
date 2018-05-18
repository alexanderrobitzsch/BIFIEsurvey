## File Name: print.object.summary.R
## File Version: 0.05

print.object.summary <- function( obji , digits ){
    obji0 <- obji
    V <- ncol(obji)
    for (vv in 1:V){
        if ( is.numeric( obji[,vv] ) ){ 
            obji[,vv] <- round( obji[,vv] , digits ) 
        }
    }
    # print(obji)
    print( obji0 , digits=digits )
}
