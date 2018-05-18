%% File Name: BIFIE.mva.Rd
%% File Version: 0.05

\name{BIFIE.mva}
\alias{BIFIE.mva}
\alias{summary.BIFIE.mva}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Missing Value Analysis
}

\description{
Conducts a missing value analysis.
}

\usage{
BIFIE.mva( BIFIEobj , missvars , covariates=NULL , se=TRUE )

\method{summary}{BIFIE.mva}(object,digits=4,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BIFIEobj}{
Object of class \code{BIFIEdata}
}
  \item{missvars}{
Vector of variables for which missing value statistics should be computed
}  
\item{covariates}{Vector of variables which work as covariates} 
\item{se}{Optional logical indicating whether statistical inference
based on replication should be employed.}
\item{object}{Object of class \code{BIFIE.correl}}
\item{digits}{Number of digits for rounding output}
\item{\dots}{Further arguments to be passed}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
A list with following entries
 \item{stat.mva}{Data frame with missing value statistics}
 \item{res_list}{List with extensive output split 
     according to each variable in \code{missvars}}
 \item{\dots}{More values}
}

%\references{
%% ~put references to the literature/web site here ~
%}

%\author{
%Alexander Robitzsch <a.robitzsch@bifie.at>
%}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%\code{\link[stats:cov.wt]{cov.wt}} (\pkg{stats}),    
%}


\examples{
#############################################################################
# EXAMPLE 1: Imputed TIMSS dataset
#############################################################################

data(data.timss1)
data(data.timssrep)

# create BIFIE.dat object
BIFIEdata <- BIFIEsurvey::BIFIE.data( data.list=data.timss1, 
                wgt=data.timss1[[1]]$TOTWGT, wgtrep=data.timssrep[, -1 ] )

# missing value analysis for "scsci" and "books" and three covariates
res1 <- BIFIEsurvey::BIFIE.mva( BIFIEdata , missvars=c("scsci" , "books" ) , 
             covariates=c("ASMMAT" , "female" , "ASSSCI") )
summary(res1)

# missing value analysis without statistical inference and without covariates
res2 <- BIFIEsurvey::BIFIE.mva( BIFIEdata, missvars=c("scsci", "books"), se=FALSE)
summary(res2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Missing value analysis}
\keyword{summary}
