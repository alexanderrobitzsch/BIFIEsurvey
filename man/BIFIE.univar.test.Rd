%% File Name: BIFIE.univar.test.Rd
%% File Version: 0.13
%% File Last Change: 2017-01-16 18:42:25

\name{BIFIE.univar.test}
\alias{BIFIE.univar.test}
\alias{summary.BIFIE.univar.test}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Analysis of Variance and Effect Sizes for Univariate Statistics
}

\description{
Computes a Wald test which tests equality of means (univariate 
analysis of variance). In addition, the \eqn{d} and \eqn{\eta}
effect sizes are computed.
}

\usage{
BIFIE.univar.test(BIFIE.method , wald_test=TRUE)

\method{summary}{BIFIE.univar.test}(object,digits=4,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BIFIE.method}{
Object of class \code{BIFIE.univar}
}
\item{wald_test}{Optional logical indicating whether a Wald test should
	 be performed.}
\item{object}{Object of class \code{BIFIE.univar.test}}
\item{digits}{Number of digits for rounding output}
\item{\dots}{Further arguments to be passed}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
A list with following entries
 \item{stat.F}{Data frame with \eqn{F} statistic for Wald test}
 \item{stat.eta}{Data frame with \eqn{\eta} effect size and its inference}
 \item{stat.dstat}{Data frame with Cohen's \eqn{d} effect size and its inference} 
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
\code{\link{BIFIE.univar}}
}


\examples{
#############################################################################
# EXAMPLE 1: Imputed TIMSS dataset - One grouping variable
#############################################################################

data(data.timss1)
data(data.timssrep)

# create BIFIE.dat object
bdat <- BIFIEsurvey::BIFIE.data( data.list=data.timss1, wgt=data.timss1[[1]]$TOTWGT,
           wgtrep=data.timssrep[, -1 ] )

#**** Model 1: 3 variables splitted by book
res1 <- BIFIEsurvey::BIFIE.univar( bdat, vars=c("ASMMAT", "ASSSCI","scsci"), 
                    group="books")
summary(res1)
# analysis of variance
tres1 <- BIFIEsurvey::BIFIE.univar.test(res1)
summary(tres1)

#**** Model 2: One variable splitted by gender
res2 <- BIFIEsurvey::BIFIE.univar( bdat , vars=c("ASMMAT") , group="female" )
summary(res2)
# analysis of variance
tres2 <- BIFIEsurvey::BIFIE.univar.test(res2)
summary(tres2)

\dontrun{
#**** Model 3: Univariate statistic: math
res3 <- BIFIEsurvey::BIFIE.univar( bdat , vars=c("ASMMAT") )
summary(res3)
tres3 <- BIFIEsurvey::BIFIE.univar.test(res3)

#############################################################################
# EXAMPLE 2: Imputed TIMSS dataset - Two grouping variables
#############################################################################

data(data.timss1)
data(data.timssrep)

# create BIFIE.dat object
bdat <- BIFIEsurvey::BIFIE.data( data.list=data.timss1 , wgt=data.timss1[[1]]$TOTWGT ,
                  wgtrep=data.timssrep[, -1 ] )

#**** Model 1: 3 variables splitted by book and female
res1 <- BIFIEsurvey::BIFIE.univar(bdat, vars=c("ASMMAT", "ASSSCI","scsci"), 
                  group=c("books","female"))
summary(res1)

# analysis of variance
tres1 <- BIFIEsurvey::BIFIE.univar.test(res1)
summary(tres1)

# extract data frame with Cohens d statistic
dstat <- tres1$stat.dstat

# extract d values for gender comparisons with same value of books
# -> 'books' refers to the first variable
ind <- which( 
  unlist( lapply( strsplit( dstat$groupval1, "#"), FUN=function(vv){vv[1]}) ) ==
  unlist( lapply( strsplit( dstat$groupval2, "#"), FUN=function(vv){vv[1]}) )
        )
dstat[ ind , ]
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Descriptive statistics}
\keyword{summary}
