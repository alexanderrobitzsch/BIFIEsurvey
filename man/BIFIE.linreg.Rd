%% File Name: BIFIE.linreg.Rd
%% File Version: 0.22
%% File Last Change: 2017-01-16 19:36:13

\name{BIFIE.linreg}
\alias{BIFIE.linreg}
\alias{summary.BIFIE.linreg}
\alias{coef.BIFIE.linreg}
\alias{vcov.BIFIE.linreg}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Linear Regression
}

\description{
Computes linear regression.
}

\usage{
BIFIE.linreg(BIFIEobj, dep=NULL , pre=NULL, formula=NULL , 
    group = NULL, group_values = NULL, se=TRUE)

\method{summary}{BIFIE.linreg}(object,digits=4,...)

\method{coef}{BIFIE.linreg}(object,...)

\method{vcov}{BIFIE.linreg}(object,...)
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
\item{object}{Object of class \code{BIFIE.linreg}}
\item{digits}{Number of digits for rounding output}
\item{\dots}{Further arguments to be passed}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
A list with following entries
 \item{stat}{Data frame with unstandardized and standardized regression
 coefficients, residual standard deviation and \eqn{R^2}}
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
\code{\link[survey:svyglm]{survey::svyglm}},
\code{\link[intsvy:timss.reg]{intsvy::timss.reg}},	
\code{\link[intsvy:timss.reg.pv]{intsvy::timss.reg.pv}},
\code{\link[stats:lm]{stats::lm}}
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

#**** Model 1: Linear regression for mathematics score
mod1 <- BIFIEsurvey::BIFIE.linreg( bdat, dep= "ASMMAT", pre=c("one","books","migrant"),
              group= "female" )
summary(mod1)

\dontrun{
# same model but specified with R formulas 
mod1a <- BIFIEsurvey::BIFIE.linreg( bdat , formula = ASMMAT ~ books + migrant ,
               group= "female" , group_values = 0:1 )
summary(mod1a)

# compare result with lm function and first imputed dataset
dat1 <- data.timss1[[1]]
mod1b <- stats::lm( ASMMAT ~ 0 + as.factor(female) + as.factor(female):books +
                              as.factor(female):migrant , 
                         data= dat1 ,  weights=dat1$TOTWGT )
summary(mod1b)

#**** Model 2: Like Model 1, but books is now treated as a factor 
mod2 <- BIFIEsurvey::BIFIE.linreg( bdat , formula = ASMMAT ~ as.factor(books) + migrant)
summary(mod2)

#############################################################################
# EXAMPLE 2: PISA data | Nonlinear regression models
#############################################################################

data(data.pisaNLD)
data <- data.pisaNLD

#--- Create BIFIEdata object immediately using BIFIE.data.jack function
bdat <- BIFIEsurvey::BIFIE.data.jack( data.pisaNLD, jktype = "RW_PISA", cdata=TRUE)
summary(bdat)

#****************************************************
#*** Model 1: linear regression
mod1 <- BIFIEsurvey::BIFIE.linreg( bdat , formula = MATH ~ HISEI )
summary(mod1)

#****************************************************
#*** Model 2: Cubic regression
mod2 <- BIFIEsurvey::BIFIE.linreg( bdat , 
              formula = MATH ~ HISEI + I(HISEI^2) + I(HISEI^3) )
summary(mod2)

#****************************************************
#*** Model 3: B-spline regression

# test with design of HISEI values
dfr <- data.frame("HISEI" = 16:90 )
des <- stats::model.frame( ~ splines::bs( HISEI , df=5 ) , dfr )
des <- des$splines
plot( dfr$HISEI , des[,1] , type="l" , pch=1 , lwd=2 , ylim=c(0,1) )
for (vv in 2:ncol(des) ){
    lines( dfr$HISEI , des[,vv] , lty=vv , col=vv , lwd=2)
}

# apply B-spline regression in BIFIEsurvey::BIFIE.linreg
mod3 <- BIFIEsurvey::BIFIE.linreg( bdat , formula = MATH ~ splines::bs(HISEI,df=5) )
summary(mod3)

#*** include transformed HISEI values for B-spline matrix in bdat
bdat2 <- BIFIEsurvey::BIFIE.data.transform( bdat , ~ 0 + splines::bs( HISEI , df = 5 ))
bdat2$varnames[ bdat2$varsindex.added ] <- paste0("HISEI_bsdes" , 
            seq( 1 , length( bdat2$varsindex.added ) ) )

#****************************************************
#*** Model 4: Nonparametric regression using BIFIE.by

?BIFIE.by

#---- (1) test function with one dataset
dat1 <- bdat$dat1
vars <- c("MATH" , "HISEI")
X <- dat1[,vars]
w <- bdat$wgt
X <- as.data.frame(X)
# estimate model
mod <- stats::loess( MATH ~ HISEI , weights = w , data = X )

# predict HISEI values
hisei_val <- data.frame( "HISEI" = seq(16,90) )
y_pred <- stats::predict( mod , hisei_val )
graphics::plot( hisei_val$HISEI , y_pred , type="l")

#--- (2) define loess function
loess_fct <- function(X,w){
    X1 <- data.frame( X , w )
    colnames(X1) <- c( vars , "wgt")
    X1 <- stats::na.omit(X1)  
#    mod <- stats::lm( MATH ~ HISEI , weights = X1$wgt , data = X1 )
    mod <- stats::loess( MATH ~ HISEI , weights = X1$wgt , data = X1 )
    y_pred <- stats::predict( mod , hisei_val )
    return(y_pred)
}

#--- (3) estimate model
mod4 <-  BIFIEsurvey::BIFIE.by( bdat , vars  , userfct= loess_fct )
summary(mod4)

# plot linear function pointwise and confidence intervals
graphics::plot( hisei_val$HISEI , mod4$stat$est , type="l" , lwd = 2 , 
        xlab="HISEI" , ylab="PVMATH" , ylim=c(430,670) ) 
graphics::lines( hisei_val$HISEI , mod4$stat$est - 1.96* mod4$stat$SE , lty=3 )
graphics::lines( hisei_val$HISEI , mod4$stat$est + 1.96* mod4$stat$SE , lty=3 )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Linear regression}
\keyword{summary}
\keyword{coef}
\keyword{vcov}

