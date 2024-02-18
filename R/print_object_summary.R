## File Name: print_object_summary.R
## File Version: 0.161

print_object_summary <- function( obji, digits)
{
    V <- ncol(obji)
    for (vv in 1L:V){
        if ( is.numeric( obji[,vv] ) ){
            obji[,vv] <- round( obji[,vv], digits=digits )
        }
    }
    print( format(obji, scientific=FALSE), digits=digits )
}
