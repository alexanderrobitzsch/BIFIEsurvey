## File Name: zzz.R
## File Version: 0.07
#  zzz.R
#
# This function is simply copied from mice package.


# on attach
.onAttach <- function(libname,pkgname){
  d <- utils::packageDescription("BIFIEsurvey")
  packageStartupMessage("|---------------------------------------------------------",
           "--------\n"  ,
        paste("| " ,d$Package," " , d$Version," (",d$Date,")",sep="") ,
        "                                       " , 
#        "\n| Maintainer: Alexander Robitzsch <robitzsch@ipn.uni-kiel.de>      " ,
        "\n| http://www.bifie.at                                             ",
        "\n|---------------------------------------------------" ,
        "--------------\n" )
    }
version <- function(pkg="BIFIEsurvey"){
  lib <- dirname( system.file(package = pkg) )
  d <- utils::packageDescription(pkg)
  return( paste(d$Package,d$Version,d$Date,lib) )
}
