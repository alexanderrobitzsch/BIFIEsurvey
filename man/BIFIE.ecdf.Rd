\name{BIFIE.ecdf}
\alias{BIFIE.ecdf}
\alias{summary.BIFIE.ecdf}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Empirical Distribution Function and Quantiles
}

\description{
Computes an empirical distribution function (and quantiles). 
If only some quantiles should
be calculated, then an appropriate vector of \code{breaks} (which are quantiles)
must be specified.
Statistical inference is not conducted for this method. 
}

\usage{
BIFIE.ecdf( BIFIEobj, vars , breaks=NULL, quanttype=1, group=NULL , group_values=NULL )

\method{summary}{BIFIE.ecdf}(object,digits=4,...)
}

% BIFIE.ecdf <- function( BIFIEobj , vars , breaks=NULL , quanttype=1 , 
%	group=NULL , group_values=NULL  ){

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
\item{quanttype}{
Type of calculation for quantiles. In case of \code{quanttype=1},
a linear interpolation is used (which is \code{type='i/n'} in
\code{\link[Hmisc:wtd.stats]{Hmisc::wtd.quantile}}), 
while for \code{quanttype=2} no interpolation is used.
}
  \item{group}{
Optional grouping variable
}
  \item{group_values}{
Optional vector of grouping values. This can be omitted and grouping
values will be determined automatically.
}
\item{object}{Object of class \code{BIFIE.ecdf}}
\item{digits}{Number of digits for rounding output}
\item{\dots}{Further arguments to be passed}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
A list with following entries
 \item{ecdf}{Data frame with probabilities and the empirical
 	 distribution function (See Examples).
 }
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
\code{\link[Hmisc:wtd.stats]{Hmisc::wtd.ecdf}},
\code{\link[Hmisc:wtd.stats]{Hmisc::wtd.quantile}}
}


\examples{
#############################################################################
# EXAMPLE 1: Imputed TIMSS dataset
#############################################################################

data(data.timss1)
data(data.timssrep)

# create BIFIE.dat object
bifieobj <- BIFIEsurvey::BIFIE.data( data.list=data.timss1, wgt=data.timss1[[1]]$TOTWGT ,
           wgtrep=data.timssrep[, -1 ] )

# ecdf
vars <- c( "ASMMAT" , "books")
group <- "female" ; group_values <- 0:1
# quantile type 1
res1 <- BIFIEsurvey::BIFIE.ecdf( bifieobj ,  vars = vars , group=group )
summary(res1)
res2 <- BIFIEsurvey::BIFIE.ecdf( bifieobj ,  vars = vars , group=group , quanttype=2)
# plot distribution function
ecdf1 <- res1$ecdf
plot( ecdf1$ASMMAT_female0 , ecdf1$yval , type="l")
plot( res2$ecdf$ASMMAT_female0 , ecdf1$yval , type="l" , lty=2)
plot( ecdf1$books_female0 , ecdf1$yval , type="l" , col="blue")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Empirical distribution function}
\keyword{Quantiles}
\keyword{summary}
%\keyword{plot}
