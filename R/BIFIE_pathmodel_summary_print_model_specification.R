## File Name: BIFIE_pathmodel_summary_print_model_specification.R
## File Version: 0.09

BIFIE_pathmodel_summary_print_model_specification <- function(object, digits)
{
    # print specified model
    cat("Syntax for path model \n")
    cat(object$lavaan.model, "\n")

    # print fixed reliabilities
    obji <- object$reliability
    if (!is.null(obji)){
        cat("Fixed reliabilities\n\n")
        print( round(obji, digits))
        cat("\n")
    }
}
