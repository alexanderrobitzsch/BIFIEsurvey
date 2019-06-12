%% File Name: BIFIE.waldtest.Rd
%% File Version: 0.33

\name{BIFIE.waldtest}
\alias{BIFIE.waldtest}
\alias{summary.BIFIE.waldtest}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Wald Tests for BIFIE Methods
}

\description{
This function performs a Wald test for objects of classes
\code{\link{BIFIE.by}},
\code{\link{BIFIE.correl}}, \code{\link{BIFIE.crosstab}}, \code{\link{BIFIE.freq}},
\code{\link{BIFIE.linreg}}, \code{\link{BIFIE.logistreg}} and \code{\link{BIFIE.univar}}.
}

\usage{
BIFIE.waldtest(BIFIE.method, Cdes, rdes, type=NULL)

\method{summary}{BIFIE.waldtest}(object,digits=4,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BIFIE.method}{
Object of classes \code{\link{BIFIE.by}},
\code{\link{BIFIE.correl}}, \code{\link{BIFIE.crosstab}},
\code{\link{BIFIE.freq}},  \code{\link{BIFIE.linreg}},
\code{\link{BIFIE.logistreg}} or \code{\link{BIFIE.univar}}
(see \code{parnames} in the Output of these methods for saved parameters)
}
\item{Cdes}{Design matrix \eqn{C} (see Details)}
 \item{rdes}{Design vector \eqn{r} (see Details)}
 \item{type}{Only applies to \code{BIFIE.correl}. In case of \code{type="cov"}
 covariances instead of correlations are used for parameter tests.}
 \item{object}{Object of class \code{BIFIE.waldtest}}
\item{digits}{Number of digits for rounding output}
\item{\dots}{Further arguments to be passed}
}

\details{
The Wald test is conducted for a parameter vector \eqn{\bold{\theta}},
specifying the hypothesis \eqn{C \bold{\theta}=r}. Statistical inference
is performed by using the \eqn{D_1} and the \eqn{D_2} statistic
(Enders, 2010, Ch. 8).

For objects of class \code{bifie.univar}, only hypotheses with respect
to means are implemented.
}


\value{
A list with following entries
 \item{stat.D}{Data frame with \eqn{D_1} and \eqn{D_2} statistic, degrees of freedom
 and p value}
 \item{\dots}{More values}
}

\references{
Enders, C. K. (2010). \emph{Applied missing data analysis}. Guilford Press.
}

%\author{
%Alexander Robitzsch <a.robitzsch@bifie.at>
%}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[survey:regTermTest]{survey::regTermTest}},
\code{\link[survey:anova.svyglm]{survey::anova.svyglm}},
\code{car::linearHypothesis}
}


\examples{
#############################################################################
# EXAMPLE 1: Imputed TIMSS dataset
#############################################################################

data(data.timss1)
data(data.timssrep)

# create BIFIE.dat object
bdat <- BIFIEsurvey::BIFIE.data( data.list=data.timss1, wgt=data.timss1[[1]]$TOTWGT,
           wgtrep=data.timssrep[, -1 ] )

#******************
#*** Model 1: Linear regression
res1 <- BIFIEsurvey::BIFIE.linreg( bdat, dep="ASMMAT", pre=c("one","books","migrant"),
         group="female" )
summary(res1)

#*** Wald test which tests whether sigma and R^2 values are the same
res1$parnames    # parameter names
pn <- res1$parnames ; PN <- length(pn)
Cdes <- matrix(0,nrow=2, ncol=PN)
colnames(Cdes) <- pn
# equality of R^2  ( R^2(female0) - R^2(female1)=0 )
Cdes[ 1, c("R^2_NA_female_0", "R^2_NA_female_1" ) ] <- c(1,-1)
# equality of sigma ( sigma(female0) - sigma(female1)=0)
Cdes[ 2, c("sigma_NA_female_0", "sigma_NA_female_1" ) ] <- c(1,-1)
# design vector
rdes <- rep(0,2)
# perform Wald test
wmod1 <- BIFIEsurvey::BIFIE.waldtest( BIFIE.method=res1, Cdes=Cdes, rdes=rdes )
summary(wmod1)

\dontrun{
#******************
#*** Model 2: Correlations

# compute some correlations
res2a <- BIFIEsurvey::BIFIE.correl( bdat, vars=c("ASMMAT","ASSSCI","migrant","books"))
summary(res2a)

# test whether r(MAT,migr)=r(SCI,migr) and r(MAT,books)=r(SCI,books)
pn <- res2a$parnames; PN <- length(pn)
Cdes <- matrix( 0, nrow=2, ncol=PN )
colnames(Cdes) <- pn
Cdes[ 1, c("ASMMAT_migrant", "ASSSCI_migrant") ] <- c(1,-1)
Cdes[ 2, c("ASMMAT_books", "ASSSCI_books") ] <- c(1,-1)
rdes <- rep(0,2)
# perform Wald test
wres2a <- BIFIEsurvey::BIFIE.waldtest( res2a, Cdes, rdes )
summary(wres2a)

#******************
#*** Model 3: Frequencies

# Number of books splitted by gender
res3a <- BIFIEsurvey::BIFIE.freq( bdat, vars=c("books"), group="female" )
summary(res3a)

# test whether book(cat4,female0)+book(cat5,female0)=book(cat4,female1)+book(cat5,female5)
pn <- res3a$parnames
PN <- length(pn)
Cdes <- matrix( 0, nrow=1, ncol=PN )
colnames(Cdes) <- pn
Cdes[ 1, c("books_4_female_0", "books_5_female_0",
    "books_4_female_1", "books_5_female_1" ) ] <- c(1,1,-1,-1)
rdes <- c(0)
# Wald test
wres3a <- BIFIEsurvey::BIFIE.waldtest( res3a, Cdes, rdes )
summary(wres3a)

#******************
#*** Model 4: Means

# math and science score splitted by gender
res4a <- BIFIEsurvey::BIFIE.univar( bdat, vars=c("ASMMAT","ASSSCI"), group="female")
summary(res4a)

# test whether there are significant gender differences in math and science
#=> multivariate ANOVA
pn <- res4a$parnames
PN <- length(pn)
Cdes <- matrix( 0, nrow=2, ncol=PN )
colnames(Cdes) <- pn
Cdes[ 1, c("ASMMAT_female_0", "ASMMAT_female_1"  ) ] <- c(1,-1)
Cdes[ 2, c("ASSSCI_female_0", "ASSSCI_female_1"  ) ] <- c(1,-1)
rdes <- rep(0,2)
# Wald test
wres4a <- BIFIEsurvey::BIFIE.waldtest( res4a, Cdes, rdes )
summary(wres4a)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Wald test}
%% \keyword{summary}
