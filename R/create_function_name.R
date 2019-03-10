## File Name: create_function_name.R
## File Version: 0.06


create_function_name <- function(pack, fun)
{
    requireNamespace(pack)
    lav_fun_00 <- NULL
    fn <- paste0(pack, paste0(rep(":",2), collapse=""), fun)
    r_op <- paste0("lav_fun_00 <- ", fn)
    eval(parse(text=r_op), envir=parent.frame())
    return(lav_fun_00)
}
