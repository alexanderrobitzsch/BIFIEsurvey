## File Name: bifie_extend_list_length2.R
## File Version: 0.03

bifie_extend_list_length2 <- function(x)
{
    N <- length(x)
    if (N==1){
        x <- list( x[[1]], x[[1]] )
    }
    return(x)
}
