%% File Name: BIFIE.logistreg.Rd
%% File Version: 0.21

\name{BIFIE.logistreg}
\alias{BIFIE.logistreg}
\alias{summary.BIFIE.logistreg}
\alias{coef.BIFIE.logistreg}
\alias{vcov.BIFIE.logistreg}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Logistic Regression
}

\description{
Computes logistic regression.
}

\usage{
BIFIE.logistreg(BIFIEobj, dep=NULL, pre=NULL, formula=NULL,
    group=NULL, group_values=NULL, se=TRUE, eps=1E-8, maxiter=100)

\method{summary}{BIFIE.logistreg}(object,digits=4,...)

\method{coef}{BIFIE.logistreg}(object,...)

\method{vcov}{BIFIE.logistreg}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BIFIEobj}{
Object of class \code{BIFIEdata}
}
\item{dep}{String for the dependent variable in the regression model}
 \item{pre}{
Vector of predictor variables. If the intercept should be included,
then use the variable \code{one} for specifying it (see Examples).
}
\item{formula}{An \R formula object which can be applied instead of
    providing \code{dep} and \code{pre}. Note that there is
   additional computation time needed for model matrix creation.}
  \item{group}{
Optional grouping variable(s)
}
  \item{group_values}{
Optional vector of grouping values. This can be omitted and grouping
values will be determined automatically.
}
\item{se}{Optional logical indicating whether statistical inference
based on replication should be employed.}
\item{eps}{Convergence criterion for parameters}
\item{maxiter}{Maximum number of iterations}
\item{object}{Object of class \code{BIFIE.logistreg}}
\item{digits}{Number of digits for rounding output}
\item{\dots}{Further arguments to be passed}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
A list with following entries
 \item{stat}{Data frame with regression coefficients}
 \item{output}{Extensive output with all replicated statistics}
 \item{\dots}{More values}
}

%\references{
%% ~put references to the literature/web site here ~
%}

%\author{
%Alexander Robitzsch <a.robitzsch@bifie.at>
%}

%---
%\section{Function Maintainer}{
%Alexander Robitzsch <a.robitzsch@bifie.at>
%}
%---

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[survey:svyglm]{survey::svyglm}},
\code{\link[stats:glm]{stats::glm}}

For linear regressions see \code{\link{BIFIE.linreg}}.
}


\examples{
#############################################################################
# EXAMPLE 1: TIMSS dataset | Logistic regression
#############################################################################

data(data.timss2)
data(data.timssrep)

# create BIFIE.dat object
bdat <- BIFIEsurvey::BIFIE.data( data.list=data.timss2, wgt=data.timss2[[1]]$TOTWGT,
                      wgtrep=data.timssrep[, -1 ] )

#**** Model 1: Logistic regression - prediction of migrational background
res1 <- BIFIEsurvey::BIFIE.logistreg( BIFIEobj=bdat, dep="migrant",
           pre=c("one","books","lang"), group="female", se=FALSE )
summary(res1)

\dontrun{
# same model, but with formula specification and standard errors
res1a <- BIFIEsurvey::BIFIE.logistreg( BIFIEobj=bdat,
              formula=migrant ~ books + lang, group="female"  )
summary(res1a)

#############################################################################
# SIMULATED EXAMPLE 2: Comparison of stats::glm and BIFIEsurvey::BIFIE.logistreg
#############################################################################

#*** (1) simulate data
set.seed(987)
N <- 300
x1 <- stats::rnorm(N)
x2 <- stats::runif(N)
ypred <- -0.75+.2*x1 + 3*x2
y <- 1*( stats::plogis(ypred) > stats::runif(N) )
data <- data.frame( "y"=y, "x1"=x1, "x2"=x2 )

#*** (2) estimation logistic regression using glm
mod1 <- stats::glm( y ~ x1 + x2, family="binomial")

#*** (3) estimation logistic regression using BIFIEdata
# create BIFIEdata object by defining 30 Jackknife zones
bifiedata <- BIFIEsurvey::BIFIE.data.jack( data, jktype="JK_RANDOM", ngr=30 )
summary(bifiedata)
# estimate logistic regression
mod2 <- BIFIEsurvey::BIFIE.logistreg( bifiedata, formula=y ~ x1+x2 )

#*** (4) compare results
summary(mod2)    # BIFIE.logistreg
summary(mod1)   # glm
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Logistic regression}
%% \keyword{summary}
%% \keyword{coef}
%% \keyword{vcov}

