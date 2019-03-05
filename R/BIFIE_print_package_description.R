## File Name: BIFIE_print_package_description.R
## File Version: 0.01

BIFIE_print_package_description <- function(pack="BIFIEsurvey")
{
    d1 <- packageDescription(pack)
    cat( paste( d1$Package, " ", d1$Version, " (", d1$Date, ")", sep=""), "\n" )
}
