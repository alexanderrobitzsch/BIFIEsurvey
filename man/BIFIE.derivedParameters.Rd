%% File Name: BIFIE.derivedParameters.Rd
%% File Version: 0.13
%% File Last Change: 2017-01-16 18:32:49

\name{BIFIE.derivedParameters}
\alias{BIFIE.derivedParameters}
\alias{summary.BIFIE.derivedParameters}
\alias{coef.BIFIE.derivedParameters}
\alias{vcov.BIFIE.derivedParameters}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Statistical Inference for Derived Parameters
}

\description{
This function performs statistical for derived parameters for objects of classes
\code{\link{BIFIE.by}}, 
\code{\link{BIFIE.correl}}, \code{\link{BIFIE.crosstab}}, \code{\link{BIFIE.freq}}, 
\code{\link{BIFIE.linreg}}, \code{\link{BIFIE.logistreg}} and \code{\link{BIFIE.univar}}.
}
%% allow different methods

\usage{
BIFIE.derivedParameters( BIFIE.method , derived.parameters , type=NULL)

\method{summary}{BIFIE.derivedParameters}(object,digits=4,...)

\method{coef}{BIFIE.derivedParameters}(object,...)

\method{vcov}{BIFIE.derivedParameters}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BIFIE.method}{
Object of classes \code{\link{BIFIE.by}}, \code{\link{BIFIE.correl}}, 
\code{\link{BIFIE.crosstab}}, 
\code{\link{BIFIE.freq}},  \code{\link{BIFIE.linreg}},
\code{\link{BIFIE.logistreg}} or \code{\link{BIFIE.univar}}
(see \code{parnames} in the Output of these methods for saved parameters)
}
\item{derived.parameters}{List with \R formulas for derived
parameters (see Examples for specification)}
\item{type}{Only applies to \code{BIFIE.correl}. In case of \code{type="cov"}
 covariances instead of correlations are used for derived parameters.}
\item{object}{Object of class \code{BIFIE.derivedParameters}}
\item{digits}{Number of digits for rounding output}
\item{\dots}{Further arguments to be passed}

}

\details{
The distribution of derived parameters is derived by the direct
calculation using original resampled parameters.
}


\value{
A list with following entries
 \item{stat}{Data frame with statistics}
 \item{coef}{Estimates of derived parameters}
 \item{vcov}{Covariance matrix of derived parameters}
 \item{parnames}{Parameter names}
 \item{\dots}{More values}
}

%\references{
%Enders, C. K. (2010). \emph{Applied missing data analysis}. Guilford Press.
%}

%\author{
%Alexander Robitzsch <a.robitzsch@bifie.at>
%}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{BIFIE.waldtest}} for multi-parameter tests.

See \code{\link[car:deltaMethod]{car::deltaMethod}} for the delta
method assuming that the multivariate distribution of the parameters is
asymptotically normal.
}
%\code{\link[survey:regTermTest]{regTermTest}} (\pkg{survey}),
%\code{\link[survey:anova.svyglm]{anova.svyglm}} (\pkg{survey})


\examples{
#############################################################################
# EXAMPLE 1: Imputed TIMSS dataset
#            Inference for correlations and derived parameters
#############################################################################

data(data.timss1)
data(data.timssrep)

# create BIFIE.dat object
bdat <- BIFIEsurvey::BIFIE.data( data.list=data.timss1, wgt=data.timss1[[1]]$TOTWGT,
           wgtrep=data.timssrep[, -1 ] )

# compute correlations
res1 <- BIFIEsurvey::BIFIE.correl( bdat , 
            vars=c("ASSSCI" , "ASMMAT" , "books" , "migrant" )  )
summary(res1)
res1$parnames
  ##    [1] "ASSSCI_ASSSCI"   "ASSSCI_ASMMAT"   "ASSSCI_books"    "ASSSCI_migrant" 
  ##    [5] "ASMMAT_ASMMAT"   "ASMMAT_books"    "ASMMAT_migrant"  "books_books"    
  ##    [9] "books_migrant"   "migrant_migrant"

# define four derived parameters
derived.parameters <- list( 
        # squared correlation of science and mathematics
        "R2_sci_mat" = ~ 0 + I( 100* ASSSCI_ASMMAT^2  )  ,
        # partial correlation of science and mathematics controlling for books
        "parcorr_sci_mat" = ~ 0 + I( ( ASSSCI_ASMMAT - ASSSCI_books * ASMMAT_books ) /
                            sqrt(( 1 - ASSSCI_books^2 ) * ( 1-ASMMAT_books^2 ) ) ) ,
        # original correlation science and mathematics (already contained in res1)
        "cor_sci_mat" = ~ 0 + I( ASSSCI_ASMMAT ) ,
        # original correlation books and migrant
        "cor_book_migra" = ~ 0 + I( books_migrant )         
        )

# statistical inference for derived parameters
res2 <- BIFIEsurvey::BIFIE.derivedParameters( res1 , derived.parameters )
summary(res2)                      
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Derived parameters}
\keyword{summary}
\keyword{coef}
\keyword{vcov}

