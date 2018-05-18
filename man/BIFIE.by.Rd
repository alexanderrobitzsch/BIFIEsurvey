%% File Name: BIFIE.by.Rd
%% File Version: 0.18

\name{BIFIE.by}
\alias{BIFIE.by}
\alias{summary.BIFIE.by}
\alias{coef.BIFIE.by}
\alias{vcov.BIFIE.by}


%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Statistics for User Defined Functions
}

\description{
Computes statistics for user defined functions.
}

\usage{
BIFIE.by( BIFIEobj , vars , userfct , userparnames=NULL ,
     group=NULL , group_values=NULL , se=TRUE , use_Rcpp = TRUE)

\method{summary}{BIFIE.by}(object,digits=4,...)

\method{coef}{BIFIE.by}(object,...)

\method{vcov}{BIFIE.by}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BIFIEobj}{
Object of class \code{BIFIEdata}
}
  \item{vars}{
Vector of variables for which statistics should be computed
}
\item{userfct}{User defined function. This function must include
a matrix \code{X} and a weight vector \code{w} as
arguments. The value of this function must be a vector.
}
\item{userparnames}{An optional vector of parameter names for the
    value of \code{userfct}.    
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
\item{use_Rcpp}{Optional logical indicating whether the user defined
function should be evaluated in \pkg{Rcpp}.}
\item{object}{Object of class \code{BIFIE.by}}
\item{digits}{Number of digits for rounding output}
\item{\dots}{Further arguments to be passed}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
A list with following entries
 \item{stat}{Data frame with statistics defined in \code{userfct}}
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
\code{\link[survey:svyby]{survey::svyby}}
}


\examples{
#############################################################################
# EXAMPLE 1: Imputed TIMSS dataset
#############################################################################

data(data.timss1)
data(data.timssrep)

# create BIFIE.dat object
bifieobj <- BIFIEsurvey::BIFIE.data( data.list=data.timss1, wgt=data.timss1[[1]]$TOTWGT,
           wgtrep=data.timssrep[, -1 ] )

#****************************           
#*** Model 1: Weighted means (as a toy example)
userfct <- function(X,w){
        pars <- c( stats::weighted.mean( X[,1] , w ) , 
                     stats::weighted.mean(X[,2] , w )   )
        return(pars)
                        }
res1 <-  BIFIEsurvey::BIFIE.by( bifieobj , vars =c("ASMMAT" , "migrant" , "books"), 
                userfct=userfct , userparnames =c("MW_MAT", "MW_Migr") , 
                group="female" )
summary(res1)

# evaluate function in pure R implementation using the use_Rcpp argument
res1b <-  BIFIEsurvey::BIFIE.by( bifieobj , vars =c("ASMMAT" , "migrant" , "books" ) , 
                userfct=userfct , userparnames =c("MW_MAT", "MW_Migr") , 
                group="female" , use_Rcpp = FALSE )
summary(res1b)

#--- statistical inference for a derived parameter (see ?BIFIE.derivedParameters)
# define gender difference for mathematics score (divided by 100)
derived.parameters <- list(
        "gender_diff" = ~ 0 + I( ( MW_MAT_female1 - MW_MAT_female0 ) / 100 )
                            )
# inference derived parameter
res1d <- BIFIEsurvey::BIFIE.derivedParameters( res1 , 
                derived.parameters = derived.parameters )
summary(res1d)

\dontrun{
#****************************
#**** Model 2: Robust linear model

# (1) start from scratch to formulate the user function for X and w
dat1 <- bifieobj$dat1
vars <- c("ASMMAT" , "migrant" , "books" )
X <- dat1[,vars]
w <- bifieobj$wgt
library(MASS)
# ASMMAT ~ migrant + books
mod <- MASS::rlm( X[,1] ~  as.matrix( X[, -1 ] ) , weights = w )
coef(mod)
# (2) define a user function "my_rlm"
my_rlm <- function(X,w){
    mod <- MASS::rlm( X[,1] ~  as.matrix( X[, -1 ] ) , weights = w )
    return( coef(mod) )
                }
# (3) estimate model
res2 <-  BIFIEsurvey::BIFIE.by( bifieobj , vars  , userfct=my_rlm , 
                group="female", group_values=0:1)
summary(res2)
# estimate model without computing standard errors
res2a <-  BIFIEsurvey::BIFIE.by( bifieobj , vars  , userfct=my_rlm, 
                group="female" , se=FALSE)
summary(res2a)

# define a user function with formula language
my_rlm2 <- function(X,w){
    colnames(X) <- vars
    X <- as.data.frame(X)
    mod <- MASS::rlm( ASMMAT ~  migrant + books , weights = w , data=X)
    return( coef(mod) )
                }
# estimate model
res2b <-  BIFIEsurvey::BIFIE.by( bifieobj, vars , userfct=my_rlm2, 
                group="female", group_values=0:1)
summary(res2b)


#****************************
#**** Model 3: Number of unique values for variables in BIFIEdata

#*** define variables for which the number of unique values should be calculated
vars <- c( "female" , "books" ,"ASMMAT" )
#*** define a user function extracting these unqiue values
userfct <- function(X,w){
        pars <- apply( X , 2 , FUN = function(vv){
                     length( unique(vv))  } )
        # Note that weights are (of course) ignored in this function
        return(pars)
                        }
#*** extract number of unique values                        
res3 <-  BIFIEsurvey::BIFIE.by( bifieobj , vars = vars , userfct=userfct , 
              userparnames = paste0( vars, "_Nunique"),  se =FALSE )
summary(res3)
  ##   Statistical Inference for User Definition Function 
  ##               parm Ncases  Nweight    est
  ##   1 female_Nunique   4668 78332.99    2.0
  ##   2  books_Nunique   4668 78332.99    5.0
  ##   3 ASMMAT_Nunique   4668 78332.99 4613.4
# number of unique values in each of the five imputed datasets 
res3$output$parsrepM
  ##        [,1] [,2] [,3] [,4] [,5]
  ##   [1,]    2    2    2    2    2
  ##   [2,]    5    5    5    5    5
  ##   [3,] 4617 4619 4614 4609 4608
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{User defined function}
\keyword{summary}
\keyword{coef}
\keyword{vcov}

