## File Name: BIFIE.summary.R
## File Version: 0.231

BIFIE.summary <- function(object, print.time=TRUE)
{
    cat("------------------------------------------------------------\n")
    BIFIE_print_package_description(pack="BIFIEsurvey")
    if ( inherits(object,"BIFIE.lavaan.survey") ){
        BIFIE_print_package_description(pack="lavaan")
        BIFIE_print_package_description(pack="lavaan.survey")
    }

    #* function
    cat( paste0("\nFunction '", class(object) ) )
    if ( inherits(object,"BIFIE.waldtest") ){
        cat( paste0( "' for BIFIE method '", object$class.BIFIE.method ) )
    }
    cat("'" )
    cat("\n\nCall:\n", paste(deparse(object$CALL), sep="\n", collapse="\n"),
                "\n\n", sep="")
    cat( "Date of Analysis:", paste( object$time[1] ), "\n" )
    if (print.time){
        cat("Computation time:", print(object$time[2] - object$time[1] ), "\n\n")
    } else {
        cat("\n")
    }

    if ( ! object$NMI ){
        cat("Multiply imputed dataset\n\n")
    }
    if ( object$NMI ){
        cat("Nested multiply imputed dataset\n\n")
    }

    cat( "Number of persons", "=", object$N, "\n" )
    # cat( "Number of imputed datasets=", object$Nimp, "\n" )
    if ( ! object$NMI){
        cat( "Number of imputed datasets", "=", object$Nimp, "\n" )
    }
    if ( object$NMI){
        cat( "Number of imputed between-nest datasets", "=", object$Nimp_NMI[1], "\n" )
        cat( "Number of imputed within-nest datasets", "=", object$Nimp_NMI[2], "\n" )
    }

    cat( "Number of Jackknife zones per dataset",  "=", object$RR, "\n" )
    cat( "Fay factor", "=", round( object$fayfac, 5 ), "\n\n" )
}
