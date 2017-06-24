\name{BIFIE.hist}
\alias{BIFIE.hist}
\alias{summary.BIFIE.hist}
\alias{plot.BIFIE.hist}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Histogram
}

\description{
Computes a histogram with same output as in 
\code{\link[graphics:hist]{graphics::hist}}.
Statistical inference is not conducted for this method.
}

\usage{
BIFIE.hist( BIFIEobj, vars , breaks=NULL, group=NULL , group_values=NULL  )

\method{summary}{BIFIE.hist}(object,...)

\method{plot}{BIFIE.hist}(x,ask=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BIFIEobj}{
Object of class \code{BIFIEdata}
}
  \item{vars}{
Vector of variables for which statistics should be computed.
}
\item{breaks}{
Optional vector of breaks. Otherwise, it will be automatically defined.
}
  \item{group}{
Optional grouping variable(s)
}
  \item{group_values}{
Optional vector of grouping values. This can be omitted and grouping
values will be determined automatically.
}
\item{object}{Object of class \code{BIFIE.hist}}
\item{x}{Object of class \code{BIFIE.hist}}
\item{ask}{Optional logical whether it should be asked for new plots.}
\item{\dots}{Further arguments to be passed}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
A list with following entries
 \item{histobj}{List with objects of class \code{histogram}}
 \item{output}{More extensive output}
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

\seealso{
\code{\link[graphics:hist]{graphics::hist}}	
}


\examples{
#############################################################################
# EXAMPLE 1: Imputed TIMSS dataset
#############################################################################

data(data.timss1)
data(data.timssrep)

# create BIFIE.dat object
bifieobj <- BIFIEsurvey::BIFIE.data( data.list=data.timss1, wgt=data.timss1[[1]]$TOTWGT,
           wgtrep=data.timssrep[, -1 ] )

# histogram  
res1 <- BIFIEsurvey::BIFIE.hist( bifieobj , vars = "ASMMAT" , group="female" )
# plot histogram for first group (female = 0)
plot( res1$histobj$ASMMAT_female0 , col="lightblue")
# plot both histograms after each other
plot( res1 )
           
# user-defined vector of breaks
res2 <- BIFIEsurvey::BIFIE.hist( bifieobj , vars = "ASMMAT", 
              breaks=seq(0,900,10), group="female" )
plot( res2 , col="orange")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Histogram}
\keyword{summary}
\keyword{plot}
