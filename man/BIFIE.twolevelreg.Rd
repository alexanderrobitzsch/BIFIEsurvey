%% File Name: BIFIE.twolevelreg.Rd
%% File Version: 0.18

\name{BIFIE.twolevelreg}
\alias{BIFIE.twolevelreg}
\alias{summary.BIFIE.twolevelreg}
\alias{coef.BIFIE.twolevelreg}
\alias{vcov.BIFIE.twolevelreg}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Two Level Regression
}

\description{
This function computes the hierarchical two level model with
random intercepts and random slopes. The full maximum
likelihood estimation is conducted by means of an 
EM algorithm (Raudenbush & Bryk, 2002).
}


\usage{
BIFIE.twolevelreg( BIFIEobj , dep , formula.fixed , formula.random , idcluster , 
   wgtlevel2 = NULL , wgtlevel1 = NULL , group=NULL , group_values=NULL , 
   recov_constraint=NULL, se=TRUE, globconv = 1E-6 , maxiter = 1000 )    

\method{summary}{BIFIE.twolevelreg}(object,digits=4,...)

\method{coef}{BIFIE.twolevelreg}(object,...)

\method{vcov}{BIFIE.twolevelreg}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BIFIEobj}{
Object of class \code{BIFIEdata}
}
\item{dep}{String for the dependent variable in the regression model}
 \item{formula.fixed}{
An \R formula for fixed effects
}
\item{formula.random}{An \R formula for random effects}
\item{idcluster}{Cluster identifier. The cluster identifiers must be
sorted in the \code{BIFIE.data} object.}
\item{wgtlevel2}{Name of Level 2 weight variable}
\item{wgtlevel1}{Name of Level 1 weight variable. This is optional.
    If it is not provided, \code{wgtlevel} is calculated from
    the total weight and \code{wgtlevel2}.}
  \item{group}{
Optional grouping variable
}
  \item{group_values}{
Optional vector of grouping values. This can be omitted and grouping
values will be determined automatically.
}
\item{recov_constraint}{Matrix for constraints of random effects covariance
matrix. The random effects are numbered according to the order in 
the specification in \code{formula.random}. The first column in
\code{recov_constraint} contains the row index in the 
covariance matrix, the second column the column index and the third column
the value to be fixed.
}
\item{se}{Optional logical indicating whether statistical inference
based on replication should be employed. In case of \code{se=FALSE},
standard errors are computed as maximum likelihood estimates under
the assumption of random sampling of level 2 clusters.
}
\item{globconv}{Convergence criterion for maximum parameter change}
\item{maxiter}{Maximum number of iterations}
\item{object}{Object of class \code{BIFIE.twolevelreg}}
\item{digits}{Number of digits for rounding output}
\item{\dots}{Further arguments to be passed}
}

\details{
The implemented random slope model can be written as 
\deqn{y_{ij} = \bold{X}_{ij} \bold{\gamma}  +  \bold{Z}_{ij} \bold{u}_j +
\varepsilon_{ij}}
where \eqn{y_{ij}} is the dependent variable, \eqn{\bold{X}_{ij}}
includes the fixed effects predictors (specified by \code{formula.fixed})
and \eqn{\bold{Z}_{ij}} includes the random effects predictors 
(specified by \code{formula.random}). The random effects \eqn{\bold{u}_j}
follow a multivariate normal distribution.
		
The function also computes a variance decomposition of explained
variance due to fixed and random effects for the within and the
between level. This variance decomposition is conducted for the predictor
matrices \eqn{\bold{X}} and \eqn{\bold{Z}}. It is assumed that
\eqn{ \bold{X}_{ij} = \bold{X}_j^B + \bold{X}_{ij}^W}.
The different sources of variance are computed by formulas as
proposed in Snijders and Bosker (2012, Ch. 7).
}


\value{
A list with following entries
 \item{stat}{Data frame with coefficients and different sources
 of variance.}
 \item{output}{Extensive output with all replicated statistics}
 \item{\dots}{More values}
}

\references{
Raudenbush, S. W., & Bryk, A. S. (2002). 
\emph{Hierarchical linear models: Applications and data analysis methods}.
Thousand Oaks: Sage.	
	
Snijders, T. A. B., & Bosker, R. J. (2012). 
\emph{Multilevel analysis: An introduction to basic and advanced multilevel 
modeling}. Thousand Oaks: Sage.
}

%\author{
%Alexander Robitzsch <a.robitzsch@bifie.at>
%}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[lme4:lmer]{lme4::lmer}}
}


\examples{
\dontrun{
library(lme4)

#############################################################################
# EXAMPLE 1: Dataset data.bifie01 | TIMSS 2011
#############################################################################

data(data.bifie01)
dat <- data.bifie01
set.seed(987)

# create dataset with replicate weights and plausible values
bdat1 <- BIFIEsurvey::BIFIE.data.jack( data=dat, jktype="JK_TIMSS", jkzone="JKCZONE", 
            jkrep="JKCREP" , wgt="TOTWGT" , pv_vars = c("ASMMAT","ASSSCI") )

# create dataset without plausible values and ignoring weights
bdat2 <- BIFIEsurvey::BIFIE.data.jack( data=dat , jktype="JK_RANDOM", ngr=10 )
# => standard errors from ML estimation

#***********************************************
# Model 1: Random intercept model

#--- Model 1a: without weights, first plausible value
mod1a <- BIFIEsurvey::BIFIE.twolevelreg( BIFIEobj=bdat2, dep="ASMMAT01", 
                formula.fixed = ~ 1, formula.random = ~ 1 , idcluster="idschool", 
                wgtlevel2 = "one", se = FALSE )
summary(mod1a)

#--- Model 1b: estimation in lme4
mod1b <- lme4::lmer( ASMMAT01 ~ 1 + ( 1 | idschool) , data = dat , REML=FALSE)
summary(mod1b)

#--- Model 1c: Like Model 1a but for five plausible values and ML inference
mod1c <- BIFIEsurvey::BIFIE.twolevelreg( BIFIEobj=bdat1, dep="ASMMAT" , 
                formula.fixed = ~ 1, formula.random = ~ 1, idcluster="idschool", 
                wgtlevel2 = "one" ,  se = FALSE )
summary(mod1c)

#--- Model 1d: weights and sampling design and all plausible values
mod1d <- BIFIEsurvey::BIFIE.twolevelreg( BIFIEobj=bdat1 , dep="ASMMAT" , 
                formula.fixed = ~ 1, formula.random = ~ 1, idcluster="idschool", 
                wgtlevel2 = "SCHWGT" )
summary(mod1d)

#***********************************************
# Model 2: Random slope model

#--- Model 2a: without weights
mod2a <- BIFIEsurvey::BIFIE.twolevelreg( BIFIEobj=bdat2 , dep="ASMMAT01" , 
                formula.fixed = ~  female +  ASBG06A , formula.random = ~ ASBG06A , 
                idcluster="idschool" , wgtlevel2 = "one" ,  se = FALSE )
summary(mod2a)

#--- Model 2b: estimation in lme4
mod2b <- lme4::lmer( ASMMAT01 ~ female +  ASBG06A + ( 1 + ASBG06A | idschool) , 
                   data = dat , REML=FALSE)
summary(mod2b)

#--- Model 2c: weights and sampling design and all plausible values
mod2c <- BIFIEsurvey::BIFIE.twolevelreg( BIFIEobj=bdat1 , dep="ASMMAT" , 
                formula.fixed = ~  female +  ASBG06A , formula.random = ~ ASBG06A , 
                idcluster="idschool" , wgtlevel2 = "SCHWGT" , maxiter=500 , se=FALSE)
summary(mod2c)

#--- Model 2d: Uncorrelated intecepts and slopes

# constraint for zero covariance between intercept and slope
recov_constraint <- matrix( c(1,2,0) , ncol=3 )
mod2d <- BIFIEsurvey::BIFIE.twolevelreg( BIFIEobj=bdat2 , dep="ASMMAT01" , 
                formula.fixed = ~ female +  ASBG06A, formula.random = ~ ASBG06A, 
                idcluster="idschool" , wgtlevel2 = "one" ,  se = FALSE,
                recov_constraint=recov_constraint )
summary(mod2d)                            

#--- Model 2e: Fixed entries in the random effects covariance matrix

# two constraints for random effects covariance
# Cov(Int, Slo)=0  # zero slope for intercept and slope
# Var(Slo)=10      # slope variance of 10
recov_constraint <- matrix( c(1,2,0,  
                      2,2,10) , ncol=3 , byrow=TRUE)
mod2e <- BIFIEsurvey::BIFIE.twolevelreg( BIFIEobj=bdat2 , dep="ASMMAT01" , 
                formula.fixed = ~  female +  ASBG06A, formula.random = ~ ASBG06A, 
                idcluster="idschool" , wgtlevel2 = "one" ,  se = FALSE,
                recov_constraint=recov_constraint )
summary(mod2e)

#############################################################################
# SIMULATED EXAMPLE 2: Two-level regression with random slopes
#############################################################################

#--- (1) simulate data
set.seed(9876)
NC <- 100    # number of clusters
Nj <- 20     # number of persons per cluster
iccx <- .4   # intra-class correlation predictor
theta <- c( 0.7 , .3 )    # fixed effects
Tmat <- diag( c(.3, .1 ) ) # variances of random intercept and slope
sig2 <- .60    # residual variance
N <- NC*Nj
idcluster <- rep( 1:NC , each=Nj )
dat1 <- data.frame("idcluster" = idcluster )
dat1$X <- rep( stats::rnorm( NC , sd = sqrt(iccx) ) , each=Nj ) + 
                 stats::rnorm( N , sd = sqrt( 1 - iccx) )
dat1$Y <- theta[1] + rep( stats::rnorm(NC , sd = sqrt(Tmat[1,1] ) ) , each=Nj ) + 
      theta[2] + rep( stats::rnorm(NC, sd=sqrt(Tmat[2,2])), each=Nj )) * dat1$X +
      stats::rnorm(N, sd= sqrt(sig2) )

#--- (2) create design object 
bdat1 <- BIFIEsurvey::BIFIE.data.jack( data=dat1, jktype="JK_GROUP", jkzone="idcluster")
summary(bdat1)

#*** Model 1: Random slope model (ML standard errors)

#- estimation using BIFIE.twolevelreg
mod1a <- BIFIEsurvey::BIFIE.twolevelreg( BIFIEobj=bdat1, dep="Y" , 
                formula.fixed = ~ 1+X, formula.random = ~ 1+X, idcluster="idcluster", 
                wgtlevel2 = "one" ,  se = FALSE )
summary(mod1a)

#- estimation in lme4
mod1b <- lme4::lmer( Y ~ X + ( 1+X | idcluster) , data = dat1 , REML = FALSE  )
summary(mod1b)

#- using Jackknife for inference
mod1c <- BIFIEsurvey::BIFIE.twolevelreg( BIFIEobj=bdat1 , dep="Y" , 
                formula.fixed = ~ 1+X , formula.random = ~ 1+X , idcluster="idcluster" , 
                wgtlevel2 = "one" ,  se = TRUE )
summary(mod1c)

# extract coefficients
coef(mod1a)
coef(mod1c)
# covariance matrix
vcov(mod1a)
vcov(mod1c)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Multilevel regression}
\keyword{summary}
\keyword{coef}
\keyword{vcov}

