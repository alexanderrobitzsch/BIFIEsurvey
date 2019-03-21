## File Name: bifie_nmi_error_message.R
## File Version: 0.04

bifie_nmi_error_message <- function(fun, NMI)
{
    if (NMI){
        n1 <- paste0("'", fun, "' cannot currently handle nested multiply imputed datasets.\n")
        stop(n1)
    }
}
