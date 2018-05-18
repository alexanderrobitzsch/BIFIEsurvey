%% File Name: BIFIE.data.select.Rd
%% File Version: 0.07

\name{BIFIEdata.select}
\alias{BIFIEdata.select}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Selection of Variables and Imputed Datasets for Objects of Class \code{BIFIEdata} 
}

\description{
This function select variables and some (or all) imputed datasets of
an object of class \code{BIFIEdata} and saves the resulting
object also of class \code{BIFIEdata}. 
}

\usage{
BIFIEdata.select(bifieobj, varnames = NULL, impdata.index = NULL)
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
or compact way, see value \code{cdata}
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
See \code{\link{BIFIE.data}} for creating \code{BIFIEdata} objects.
}
%% ~~objects to See Also as \code{\link{help}}, ~~~

\examples{
#############################################################################
# EXAMPLE 1: Some manipulations of BIFIEdata objects created from data.timss1
#############################################################################
data(data.timss1)
data(data.timssrep)

# create BIFIEdata
bdat1 <- BIFIEsurvey::BIFIE.data( data.list=data.timss1 , wgt=data.timss1[[1]]$TOTWGT,
            wgtrep=data.timssrep[, -1 ])
summary(bdat1)           

# create BIFIEcdata object
bdat2 <- BIFIEsurvey::BIFIE.data( data.list=data.timss1 , wgt=data.timss1[[1]]$TOTWGT,
            wgtrep=data.timssrep[, -1 ] , cdata=TRUE )
summary(bdat2)           

# selection of variables for BIFIEdata object
bdat1a <- BIFIEsurvey::BIFIEdata.select( bdat1 , varnames = bdat1$varnames[ 1:7 ] )
# selection of variables and 1st, 2nd and 4th imputed datasets of BIFIEcdata object
bdat2a <- BIFIEsurvey::BIFIEdata.select( bdat2 , varnames = bdat2$varnames[ 1:7 ] , 
                impdata.index = c(1,2,4) )
summary(bdat1a)
summary(bdat2a)    
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{BIFIEdata}
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
