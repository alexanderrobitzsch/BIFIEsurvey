%% File Name: BIFIE.crosstab.Rd
%% File Version: 0.09
%% File Last Change: 2017-01-16 18:24:41

\name{BIFIE.crosstab}
\alias{BIFIE.crosstab}
\alias{summary.BIFIE.crosstab}
\alias{coef.BIFIE.crosstab}
\alias{vcov.BIFIE.crosstab}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Cross Tabulation
}

\description{
Creates cross tabulations and computes some effect sizes.
}

\usage{
BIFIE.crosstab( BIFIEobj , vars1 , vars2 , vars_values1 = NULL , vars_values2 = NULL , 
     group=NULL , group_values=NULL , se=TRUE )

\method{summary}{BIFIE.crosstab}(object,digits=3,...)

\method{coef}{BIFIE.crosstab}(object,...)

\method{vcov}{BIFIE.crosstab}(object,...)
}



%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BIFIEobj}{
Object of class \code{BIFIEdata}
}
  \item{vars1}{
Row variable
}
  \item{vars2}{
Column variable
}
\item{vars_values1}{Optional vector of values of variable \code{vars1}}
\item{vars_values2}{Optional vector of values of variable \code{vars2}}
  \item{group}{
Optional grouping variable(s)
}
  \item{group_values}{
Optional vector of grouping values. This can be omitted and grouping
values will be determined automatically.
}
\item{se}{Optional logical indicating whether statistical inference
based on replication should be employed.}
\item{object}{Object of class \code{BIFIE.univar}}
\item{digits}{Number of digits for rounding output}
\item{\dots}{Further arguments to be passed}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
A list with following entries
 \item{stat.probs}{Statistics for joint and conditional probabilities}
 \item{stat.marg}{Statistics for marginal probabilities}
 \item{stat.es}{Statistics for effect sizes \eqn{w} (based on \eqn{\chi^2}),
     Cramers \eqn{V}, Goodman's gamma, the PRE lambda measure and
     Kruskals tau.}
\item{chisq.test}{Table with chi square statistics (ignoring weights and
sampling scheme)}
 \item{output}{Extensive output with all replicated statistics}
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
\code{\link[survey:svytable]{survey::svytable}},	
\code{\link[Hmisc:wtd.stats]{Hmisc::wtd.table}}
}
%\code{\link[intsvy:timss.table]{timss.table}} (\pkg{intsvy}),
%\code{\link[survey:svytable]{svytable}} (\pkg{survey}),

\examples{
#############################################################################
# EXAMPLE 1: Imputed TIMSS dataset
#############################################################################

data(data.timss1)
data(data.timssrep)

# create BIFIE.dat object
bifieobj <- BIFIEsurvey::BIFIE.data( data.list=data.timss1 , wgt=data.timss1[[1]]$TOTWGT,
           wgtrep=data.timssrep[, -1 ] )

#--- Model 1: cross tabulation           
res1 <- BIFIEsurvey::BIFIE.crosstab( bifieobj , vars1="migrant" , 
               vars2 = "books", group="female" )
summary(res1)                  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Descriptive statistics}
\keyword{summary}
\keyword{coef}
\keyword{vcov}

