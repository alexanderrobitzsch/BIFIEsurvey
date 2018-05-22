%% File Name: BIFIE.data.Rd
%% File Version: 0.42

\name{BIFIE.data}
\alias{BIFIE.data}
\alias{summary.BIFIEdata}
\alias{print.BIFIEdata}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Creates an Object of Class \code{BIFIEdata}
}

\description{
This function creates an object of class \code{BIFIEdata}.
Finite sampling correction of statistical inferences can be
conducted by specifying appropriate input in the \code{fayfac}
argument.
}
\usage{
BIFIE.data(data.list, wgt=NULL, wgtrep=NULL, fayfac=1, pv_vars=NULL,
     pvpre=NULL, cdata=FALSE, NMI=FALSE)

\method{summary}{BIFIEdata}(object,...)

\method{print}{BIFIEdata}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.list}{
List of multiply imputed datasets. Can be also a list of list of imputed
datasets in case of nested multiple imputation. Then, the argument
\code{NMI=TRUE} must be specified.
}
  \item{wgt}{
A string indicating the label of case weight or a vector
containing all case weights.
}
  \item{wgtrep}{
Optional vector of replicate weights
}
\item{fayfac}{Fay factor for calculating standard errors, a numeric value.
    If finite sampling correction
    is requested, an appropriate vector input can be used (see Example 3).
}
\item{pv_vars}{Optional vector for names of plausible values, see
\code{\link{BIFIE.data.jack}}.}
\item{pvpre}{Optional vector for prefixes of plausible values, see
\code{\link{BIFIE.data.jack}}.}
\item{cdata}{An optional logical indicating whether the \code{BIFIEdata}
object should be compactly saved. The default is \code{FALSE}.}
\item{NMI}{Optional logical indicating whether \code{data.list} is obtained
by nested multiple imputation.}
\item{object}{Object of class \code{BIFIEdata}}
\item{x}{Object of class \code{BIFIEdata}}
\item{\dots}{Further arguments to be passed}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
An object of class \code{BIFIEdata} saved in a non-compact
or compact way, see value \code{cdata}. The following entries are
included in the list:

\item{datalistM}{Stacked list of imputed datasets (if \code{cdata=FALSE})}
\item{wgt}{Vector with case weights}
\item{wgtrep}{Matrix with replicate weights}
\item{Nimp}{Number of imputed datasets}
\item{N}{Number of observations in a dataset}
\item{dat1}{Last imputed dataset}
\item{varnames}{Vector with variable names}
\item{fayfac}{Fay factor.}
\item{RR}{Number of replicate weights}
\item{NMI}{Logical indicating whether the dataset is nested multiply
imputed.}
\item{cdata}{Logical indicating whether the \code{BIFIEdata} object is in
    compact format (\code{cdata=TRUE}) or in a non-compact format
    (\code{cdata=FALSE}).
}
\item{Nvars}{Number of variables}
\item{variables}{Data frame including some informations about variables.
All transformations are saved in the column \code{source}.}
\item{datalistM_ind}{Data frame with response indicators
(if \code{cdata=TRUE})}
\item{datalistM_imputed}{Data frame with imputed values
(if \code{cdata=TRUE})}
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

\seealso{
See \code{\link{BIFIE.data.transform}} for data transformations on
\code{BIFIEdata} objects.

For saving and loading \code{BIFIEdata} objects see
\code{\link{save.BIFIEdata}}.

For converting PIRLS/TIMSS or PISA datasets into \code{BIFIEdata}
objects see \code{\link{BIFIE.data.jack}}.

See the \code{\link{BIFIEdata2svrepdesign}} function for converting
\code{BIFIEdata} objects to objects used in the \pkg{survey} package.
}

\examples{
#############################################################################
# EXAMPLE 1: Create BIFIEdata object with multiply-imputed TIMSS data
#############################################################################
data(data.timss1)
data(data.timssrep)

bdat <- BIFIEsurvey::BIFIE.data( data.list=data.timss1, wgt=data.timss1[[1]]$TOTWGT,
            wgtrep=data.timssrep[, -1 ] )
summary(bdat)
# create BIFIEdata object in a compact way
bdat2 <- BIFIEsurvey::BIFIE.data( data.list=data.timss1, wgt=data.timss1[[1]]$TOTWGT,
            wgtrep=data.timssrep[, -1 ], cdata=TRUE)
summary(bdat2)

\dontrun{
#############################################################################
# EXAMPLE 2: Create BIFIEdata object with one dataset
#############################################################################
data(data.timss2)

# use first dataset with missing data from data.timss2
bdat <- BIFIEsurvey::BIFIE.data( data.list=data.timss2[[1]], wgt=data.timss2[[1]]$TOTWGT)
}

#############################################################################
# EXAMPLE 3: BIFIEdata objects with finite sampling correction
#############################################################################

data(data.timss1)
data(data.timssrep)

#-----
# BIFIEdata object without finite sampling correction
bdat1 <- BIFIEsurvey::BIFIE.data( data.list=data.timss1, wgt=data.timss1[[1]]$TOTWGT,
            wgtrep=data.timssrep[, -1 ] )
summary(bdat1)

#-----
# generate BIFIEdata object with finite sampling correction by adjusting
# the "fayfac" factor
bdat2 <- bdat1
#-- modify "fayfac" constant
fayfac0 <- bdat1$fayfac
# set fayfac=.75 for the first 50 replication zones (25% of students in the
# population were sampled) and fayfac=.20 for replication zones 51-75
# (meaning that 80% of students were sampled)
fayfac <- rep( fayfac0, bdat1$RR )
fayfac[1:50] <- fayfac0 * .75
fayfac[51:75] <- fayfac0 * .20
# include this modified "fayfac" factor in bdat2
bdat2$fayfac <- fayfac
summary(bdat2)
summary(bdat1)

#---- compare some univariate statistics
# no finite sampling correction
res1 <- BIFIEsurvey::BIFIE.univar( bdat1, vars="ASMMAT")
summary(res1)
# finite sampling correction
res2 <- BIFIEsurvey::BIFIE.univar( bdat2, vars="ASMMAT")
summary(res2)

\dontrun{
#############################################################################
# EXAMPLE 4: Create BIFIEdata object with nested multiply imputed dataset
#############################################################################

data(data.timss4)
data(data.timssrep)

# nested imputed dataset, save it in compact format
bdat <- BIFIEsurvey::BIFIE.data( data.list=data.timss4,
            wgt=data.timss4[[1]][[1]]$TOTWGT, wgtrep=data.timssrep[, -1 ],
            NMI=TRUE, cdata=TRUE )
summary(bdat)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{BIFIEdata}
\keyword{summary}
\keyword{print}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
