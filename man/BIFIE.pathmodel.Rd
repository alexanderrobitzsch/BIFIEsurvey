%% File Name: BIFIE.pathmodel.Rd
%% File Version: 0.15

\name{BIFIE.pathmodel}
\alias{BIFIE.pathmodel}
\alias{summary.BIFIE.pathmodel}
\alias{coef.BIFIE.pathmodel}
\alias{vcov.BIFIE.pathmodel}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Path Model Estimation
}

\description{
This function computes a path model. Predictors are allowed to possess 
measurement errors.
Known measurement error variances (and covariances) or reliabilities
can be specified by the user. Alternatively, a set of indicators can be
defined for each latent variable, and for each imputed and replicated
dataset the measurement error variance is determined by means
of calculating the reliability Cronbachs alpha. Measurement errors
are handled by adjusting covariance matrices (see Buonaccorsi, 2010, Ch. 5).
}

\usage{
BIFIE.pathmodel( BIFIEobj , lavaan.model , reliability=NULL , group=NULL , 
        group_values=NULL , se=TRUE )
        
\method{summary}{BIFIE.pathmodel}(object,digits=4,...)

\method{coef}{BIFIE.pathmodel}(object,...)

\method{vcov}{BIFIE.pathmodel}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BIFIEobj}{
Object of class \code{BIFIEdata}
}
\item{lavaan.model}{String including the model specification in
   \pkg{lavaan} syntax. \code{lavaan.model} also allows the
   extended functionality in the 
   \code{\link[TAM:lavaanify.IRT]{TAM::lavaanify.IRT}} function.
        }
 \item{reliability}{
Optional vector containing the reliabilities of each variable.
This vector can also include only a subset of all variables.
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
\item{object}{Object of class \code{BIFIE.pathmodel}}
\item{digits}{Number of digits for rounding output}
\item{\dots}{Further arguments to be passed}
}

\details{
The following conventions are used as parameter labels in the output.

\code{Y~X} is the regression coefficient of the regression from \eqn{Y}
on \eqn{X}. 

\code{X->Z->Y} denotes the path coefficient from \eqn{X} to \eqn{Y}
passing the mediating variable \eqn{Z}.

\code{X-+>Y} denotes the total effect (of all paths) from \eqn{X} to \eqn{Y}.

\code{X-~>Y} denotes the sum of all indirect effects from \eqn{X} to \eqn{Y}.

The parameter suffix \code{_stand} refers to parameters for which
all variables are standardized.
}


\value{
A list with following entries
 \item{stat}{Data frame with unstandardized and standardized regression
 coefficients, path coefficients, total and indirect effects, 
 residual variances, and \eqn{R^2}}
 \item{output}{Extensive output with all replicated statistics}
 \item{\dots}{More values}
}

\references{
Buonaccorsi, J. P. (2010). 
\emph{Measurement error: Models, methods, and applications}. CRC Press.
}

%\author{
%Alexander Robitzsch <a.robitzsch@bifie.at>
%}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See the \pkg{lavaan} and \pkg{lavaan.survey} package.

For the \code{lavaan} syntax, see
\code{\link[lavaan:lavaanify]{lavaan::lavaanify}} and
\code{\link[TAM:lavaanify.IRT]{TAM::lavaanify.IRT}}
}


\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Path model data.bifie01
#############################################################################    

data(data.bifie01)
dat <- data.bifie01
# create dataset with replicate weights and plausible values
bifieobj <- BIFIEsurvey::BIFIE.data.jack( data=dat,  jktype="JK_TIMSS", 
                jkzone = "JKCZONE", jkrep="JKCREP", wgt="TOTWGT", 
                pv_vars = c("ASMMAT","ASSSCI") )

#**************************************************************            
#*** Model 1: Path model
lavmodel1 <- "
     ASMMAT ~ ASBG07A + ASBG07B  + ASBM03 + ASBM02A + ASBM02E
     # define latent variable with 2nd and 3rd item in reversed scoring
     ASBM03 =~ 1*ASBM03A + (-1)*ASBM03B + (-1)*ASBM03C + 1*ASBM03D
     ASBG07A ~ ASBM02E
     ASBG07A ~~ .2*ASBG07A    # measurement error variance of .20
     ASBM02E ~~ .45*ASBM02E     # measurement error variance of .45 
     ASBM02E ~ ASBM02A + ASBM02B
        "
#--- Model 1a: model calculated by gender
mod1a <- BIFIEsurvey::BIFIE.pathmodel( bifieobj , lavmodel1 , group = "female" )
summary(mod1a)
        
#--- Model 1b: Input of some known reliabilities
reliability <- c( "ASBM02B" = .6 , "ASBM02A" = .8 )
mod1b <- BIFIEsurvey::BIFIE.pathmodel( bifieobj, lavmodel1, reliability=reliability)
summary(mod1b)        
        
#**************************************************************            
#*** Model 2: Linear regression with errors in predictors

# specify lavaan model
lavmodel2 <- "
     ASMMAT ~ ASBG07A + ASBG07B + ASBM03A 
     ASBG07A ~~ .2*ASBG07A
        "
mod2 <- BIFIEsurvey::BIFIE.pathmodel( bifieobj, lavmodel2  )
summary(mod2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Linear regression}
\keyword{summary}
\keyword{coef}
\keyword{vcov}

