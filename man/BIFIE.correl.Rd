%% File Name: BIFIE.correl.Rd
%% File Version: 0.16
%% File Last Change: 2017-01-16 18:22:10

\name{BIFIE.correl}
\alias{BIFIE.correl}
\alias{summary.BIFIE.correl}
\alias{coef.BIFIE.correl}
\alias{vcov.BIFIE.correl}


%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Correlations and Covariances
}

\description{
Computes correlations and covariances
}

\usage{
BIFIE.correl(BIFIEobj, vars, group = NULL, group_values = NULL, se=TRUE)

\method{summary}{BIFIE.correl}(object,digits=4, ...)

\method{coef}{BIFIE.correl}(object,type=NULL, ...)

\method{vcov}{BIFIE.correl}(object,type=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BIFIEobj}{
Object of class \code{BIFIEdata}
}
  \item{vars}{
Vector of variables for which statistics should be computed
}
  \item{group}{
Optional grouping variable(s)
}
  \item{group_values}{
Optional vector of grouping values. This can be omitted and grouping
values will be determined automatically.
}
\item{se}{Optional logical indicating whether statistical inference
based on replication should be employed.}
\item{object}{Object of class \code{BIFIE.correl}}
\item{digits}{Number of digits for rounding output}
\item{type}{If \code{type="cov"}, then covariances instead of
correlations are extracted.}
\item{\dots}{Further arguments to be passed}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
A list with following entries
 \item{stat.cor}{Data frame with correlation statistics}
 \item{stat.cov}{Data frame with covariance statistics}
 \item{cor_matrix}{List of estimated correlation matrices}
 \item{cov_matrix}{List of estimated covariance matrices} 
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
\code{\link[stats:cov.wt]{stats::cov.wt}},	
\code{\link[intsvy:timss.rho]{intsvy::timss.rho}},
\code{\link[intsvy:timss.rho.pv]{intsvy::timss.rho.pv}},
\code{\link[Hmisc:rcorr]{Hmisc::rcorr}},
\code{\link[miceadds:ma.wtd.statNA]{miceadds::ma.wtd.corNA}}
}


\examples{
#############################################################################
# EXAMPLE 1: Imputed TIMSS dataset
#############################################################################

data(data.timss1)
data(data.timssrep)

# create BIFIE.dat object
bdat <- BIFIEsurvey::BIFIE.data( data.list=data.timss1 , wgt=  data.timss1[[1]]$TOTWGT ,
           wgtrep=data.timssrep[, -1 ] )

# Correlations splitted by gender
res1 <- BIFIEsurvey::BIFIE.correl( bdat , vars=c("lang" , "books" , "migrant" ) , 
              group="female" , group_values=0:1 )
summary(res1)  

# Correlations splitted by gender: no statistical inference (se=FALSE)
res1a <- BIFIEsurvey::BIFIE.correl( bdat , vars=c("lang" , "books" , "migrant" ) , 
              group="female" , group_values=0:1 , se=FALSE)
summary(res1a)  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Descriptive statistics}
\keyword{summary}
\keyword{coef}
\keyword{vcov}
