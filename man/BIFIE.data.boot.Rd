%% File Name: BIFIE.data.boot.Rd
%% File Version: 0.13

\name{BIFIE.data.boot}
\alias{BIFIE.data.boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create \code{BIFIE.data} Object based on Bootstrap
}

\description{
Creates a \code{BIFIE.data} object based on bootstrap designs.
The sampling is done assuming independence of cases.
}

\usage{
BIFIE.data.boot( data, wgt=NULL,  pv_vars=NULL,
         Nboot=500, seed=.Random.seed, cdata=FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data frame: Can be a single or a list of multiply imputed datasets
}
  \item{wgt}{
A string indicating the label of case weight.
}
  \item{pv_vars}{
An optional vector of plausible values which define multiply
imputed datasets.
}
\item{Nboot}{Number of
bootstrap samples for usage}
\item{seed}{Simulation seed.}
\item{cdata}{An optional logical indicating whether the \code{BIFIEdata}
object should be compactly saved. The default is \code{FALSE}.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
Object of class \code{BIFIEdata}
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
\code{\link{BIFIE.data}}, \code{\link{BIFIE.data.jack}}
}



\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Bootstrap TIMSS data set
#############################################################################
data(data.timss1)

# bootstrap samples using weights
bifieobj1 <- BIFIEsurvey::BIFIE.data.boot( data.timss1, wgt="TOTWGT" )
summary(bifieobj1)

# bootstrap samples without weights
bifieobj2 <- BIFIEsurvey::BIFIE.data.boot( data.timss1  )
summary(bifieobj2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{BIFIEdata}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
