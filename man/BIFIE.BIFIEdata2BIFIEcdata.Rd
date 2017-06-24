\name{BIFIE.BIFIEdata2BIFIEcdata}
\alias{BIFIE.BIFIEdata2BIFIEcdata}
\alias{BIFIE.BIFIEcdata2BIFIEdata}
\alias{BIFIE.BIFIEdata2datalist}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Conversion and Selection of \code{BIFIEdata} Objects
}

\description{
Functions for converting and selecting objects of class \code{BIFIEdata}.
The function \code{BIFIE.BIFIEdata2BIFIEcdata} converts the \code{BIFIEdata}
objects in a non-compact form (\code{cdata=FALSE}) into an object of
class \code{BIFIEdata} in a compact form (\code{cdata=TRUE}).
The function \code{BIFIE.BIFIE2data2BIFIEdata} takes the reverse operation.

The function \code{BIFIE.BIFIEdata2datalist} converts a (part) of the
object of class \code{BIFIEdata} into a list of multiply-imputed 
datasets.
}

\usage{
BIFIE.BIFIEdata2BIFIEcdata(bifieobj, varnames = NULL, impdata.index = NULL)

BIFIE.BIFIEcdata2BIFIEdata(bifieobj, varnames = NULL, impdata.index = NULL)

BIFIE.BIFIEdata2datalist(bifieobj, varnames = NULL, impdata.index = NULL)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bifieobj}{
Object of class \code{BIFIEdata}
}
  \item{varnames}{
Variables chosen for the selection
}
  \item{impdata.index}{
Selected indices of imputed datasets
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
An object of class \code{BIFIEdata} saved in a non-compact
or compact way, see value \code{cdata}.
}


%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{BIFIE.data}}
}


\examples{
#############################################################################
# EXAMPLE 1: BIFIEdata conversions using data.timss1 dataset
#############################################################################
data(data.timss1)
data(data.timssrep)

# create BIFIEdata object
bdat1 <- BIFIEsurvey::BIFIE.data( data.list=data.timss1, wgt= data.timss1[[1]]$TOTWGT,
            wgtrep=data.timssrep[, -1 ])
summary(bdat1)           

# convert BIFIEdata object bdat1 into a BIFIEcdata object with
#  only using the first three datasets and a variable selection
bdat2 <- BIFIEsurvey::BIFIE.BIFIEdata2BIFIEcdata( bifieobj=bdat1, 
                varnames=bdat1$varnames[ c(1:7,10) ] )

# convert bdat2 into BIFIEdata object and only use the first three imputed datasets
bdat3 <- BIFIEsurvey::BIFIE.BIFIEcdata2BIFIEdata( bifieobj=bdat2 , impdata.index=1:3)

# object summaries
summary(bdat1)
summary(bdat2)
summary(bdat3)		

\dontrun{
#############################################################################
# EXAMPLE 2: Extract unique elements in BIFIEdata object
#############################################################################

data(data.timss1)
data(data.timssrep)

# create BIFIEdata object
bifieobj <- BIFIEsurvey::BIFIE.data( data.list=data.timss1, wgt= data.timss1[[1]]$TOTWGT,
            wgtrep=data.timssrep[, -1 ])
summary(bifieobj)

# define variables for which unique values should be extracted
vars <- c( "female" , "books" ,"ASMMAT" )
# convert these variables from BIFIEdata object into a list of datasets
bdatlist <- BIFIEsurvey::BIFIE.BIFIEdata2datalist( bifieobj , varnames = vars )
# look for unique values in first dataset for variables
values <- lapply( bdatlist[[1]] , FUN = function(vv){
                sort( unique( vv ) ) } )
# number of unique values in first dataset
Nvalues <- lapply( bdatlist[[1]] , FUN = function(vv){
                length( unique( vv ) ) } )
# number of unique values in all datasets
Nvalues2 <- lapply( vars , FUN = function(vv){
    #vv <- vars[1]
    unlist( lapply( bdatlist , FUN = function(dd){
                length( unique( dd[,vv]  ) )
                        }    )     )
                    } )
# --> for extracting the number of unique values using BIFIE.by and a user
#     defined function see Example 1, Model 3 in "BIFIE.by" 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{BIFIEdata}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
