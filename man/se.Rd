%% File Name: se.Rd
%% File Version: 0.05

\name{se}
\alias{se}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Standard Errors of Estimated Parameters
}

\description{
Outputs vector of standard errors of an estimated parameter vector.
}

\usage{
se(object)
}

\arguments{
  \item{object}{
Object for which S3 method \code{vcov} can be applied
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Vector
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
\code{\link[survey:SE]{survey::SE}}
}

\examples{
#############################################################################
# EXAMPLE 1: Toy example with lm function
#############################################################################

set.seed(906)
N <- 100
x <- seq(0,1,length=N)
y <- .6*x + stats::rnorm(N, sd=1)
mod <- stats::lm( y ~ x )
coef(mod)
vcov(mod)
se(mod)
summary(mod)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{table}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
