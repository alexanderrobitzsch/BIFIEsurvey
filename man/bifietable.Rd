%% File Name: bifietable.Rd
%% File Version: 1.05

\name{bifietable}
\alias{bifietable}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
An \pkg{Rcpp} Based Version of the \code{table} Function
}

\description{
This is an \pkg{Rcpp} based version of the 
\code{\link[base:table]{base::table}} function.
}

\usage{
bifietable(vec, sort.names=FALSE)
}

\arguments{
  \item{vec}{
A numeric or character vector
}
\item{sort.names}{An optional logical indicating whether values in the
character vector should also be sorted in the table output}    
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Same output like \code{\link[base:table]{base::table}}
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
\code{\link[base:table]{base::table}}
}

\examples{
data(data.timss1)
table( data.timss1[[1]][,"books"] )
BIFIEsurvey::bifietable( data.timss1[[1]][,"books"] )    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{table}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
