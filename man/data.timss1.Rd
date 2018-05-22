%% File Name: data.timss1.Rd
%% File Version: 0.30

\name{data.timss}
\alias{data.timss}
\alias{data.timss1}
\alias{data.timss2}
\alias{data.timssrep}
\alias{data.timss3}
\alias{data.timss1.ind}
\alias{data.timss4}

\docType{data}
\title{
Dataset TIMSS 2011
}

\description{
Example dataset TIMSS 2011
}

\usage{
data(data.timss1)
data(data.timss1.ind)
data(data.timss2)
data(data.timssrep)
data(data.timss3)
data(data.timss4)
}


\format{
The dataset \code{data.timss1} is a list containing 5 imputed datasets.
The dataset \code{data.timss1.ind} contains response indicators of these
5 imputed datasets in \code{data.timss1}.

The dataset \code{data.timss2} is a list containing 5 datasets in which
only plausible values are imputed, but student covariates are missing.

The dataset \code{data.timssrep} contains replicate weights of students.

The dataset \code{data.timss3} is a TIMSS dataset with some missing
student covariates and all 5 plausible values contained in one file.

The dataset \code{data.timss4} is a list containing nested multiply imputed
datasets, with 5 between-nest and 4 within-nest imputations.
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
\dontrun{
library(survey)
library(lavaan.survey)
library(intsvy)
library(mitools)

#############################################################################
# EXAMPLE 1: TIMSS dataset data.timss3 (one dataset including all PVs)
#############################################################################

data(data.timss2)
data(data.timss3)
data(data.timssrep)

# Analysis based on official 'single' datasets (data.timss3)
# There are 5 plausible values, but student covariates are not imputed.

#--- create object of class BIFIE data
bdat3 <- BIFIEsurvey::BIFIE.data(data.timss3, wgt=data.timss3$TOTWGT,
              wgtrep=data.timssrep[,-1], fayfac=1)
summary(bdat3)
# This BIFIEdata object contains one dataset in which all
# plausible values are included. This object can be used
# in analysis without plausible values.
# Equivalently, one can define bdat3 much simpler by
bdat3 <- BIFIEsurvey::BIFIE.data.jack(data.timss3, jktype="JK_TIMSS")
summary(bdat3)

#--- In the following, the object bdat4 is defined with 5 datasets
# referring to 5 plausible values.
bdat4 <- BIFIEsurvey::BIFIE.data.jack(data.timss3, pv_vars=c("ASMMAT","ASSSCI"),
               jktype="JK_TIMSS")
summary(bdat4)

#--- create object in survey package
dat3a <- as.data.frame( cbind( data.timss2[[1]], data.timssrep ) )
RR <- ncol(data.timssrep) - 1       # number of jackknife zones
svydes3 <- survey::svrepdesign(data=dat3a, weights=~TOTWGT, type="JKn",
                 repweights='w_fstr[0-9]', scale=1,  rscales=rep(1,RR), mse=TRUE)
summary(svydes3)

#--- create object with imputed datasets in survey
datL <- data.timss2
# include replicate weights in each dataset
for (ii in 1:5){
    dat1 <- datL[[ii]]
    dat1 <- cbind(  dat1, data.timssrep[,-1] )
    datL[[ii]] <- dat1
}
datL <- mitools::imputationList(list( datL[[1]],datL[[2]],datL[[3]],datL[[4]],datL[[5]]))
svydes4 <- survey::svrepdesign(data=datL, weights=~TOTWGT, type="JKn",
                   repweights='w_fstr[0-9]', scale=1,  rscales=rep(1,RR), mse=TRUE)
summary(svydes4)

#--- reconstruct data.timss3 for intsvy package. Plausible values must be labeled
# as PV01, PV02, ... and NOT PV1, PV2, ...
data.timss3a <- data.timss3
colnames(data.timss3a) <- gsub( "ASMMAT", "ASMMAT0", colnames(data.timss3a) )
colnames(data.timss3a) <- gsub( "ASSSCI", "ASSSCI0", colnames(data.timss3a) )

#***************************
# Model 1: Linear regression (no grouping variable)

#--- linear regression in survey
mod1a <-  survey::svyglm( scsci ~ migrant + books, design=svydes3)
summary(mod1a)

#--- regression with pirls.reg (intsvy)
mod1b <- intsvy::pirls.reg( y="scsci", x=c("migrant", "books" ), data=data.timss3)
mod1b

#---- regression with BIFIEsurvey
mod1c <- BIFIEsurvey::BIFIE.linreg( bdat3, dep="scsci", pre=c("one","migrant","books"))
summary(mod1c)

#--- regression with lavaan.survey package
lavmodel <- "
    scsci ~ migrant + books
    scsci ~ 1
    scsci ~~ scsci
        "
# fit in lavaan
lavaan.fit <- lavaan::lavaan( lavmodel, data=data.timss3, estimator="MLM")
summary(lavaan.fit)
# using all replicated weights
mod1d <- lavaan.survey::lavaan.survey(lavaan.fit=lavaan.fit, survey.design=svydes3 )
summary(mod1d)

#***************************
# Model 2: Linear regression (grouped by female)

#--- linear regression in survey
mod2a <- survey::svyglm( scsci ~ 0 + as.factor(female) + as.factor(female):migrant
               + as.factor(female):books, design=svydes3)
summary(mod2a)

#--- regression with pirls.reg (intsvy)
mod2b <- intsvy::pirls.reg( y="scsci", x=c("migrant", "books" ),
                 by="female", data=data.timss3)
mod2b[["0"]]   # regression coefficients female=0
mod2b[["1"]]   # regression coefficients female=1

#--- regression with BIFIEsurvey
mod2c <- BIFIEsurvey::BIFIE.linreg( bdat3, dep="scsci",
               pre=c("one","migrant","books"), group="female")
summary(mod2c)

#--- regression with lavaan.survey package
lavmodel <- "
    scsci ~ migrant + books
    scsci ~ 1
    scsci ~~ scsci
        "
# fit in lavaan
lavaan.fit <- lavaan::lavaan( lavmodel, data=data.timss3, group="female", estimator="MLM")
summary(lavaan.fit)
mod2d <- lavaan.survey::lavaan.survey(lavaan.fit=lavaan.fit, survey.design=svydes3 )
summary(mod2d)

#***************************
# Model 3: Linear regression with mathematics PVs
library(mitools)

#--- linear regression in survey
mod3a <- with(svydes4,  survey::svyglm( ASMMAT ~ migrant + books, design=svydes4 ) )
res3a <- mitools::MIcombine(mod3a)
summary(res3a)

#--- regression with pirls.reg.pv (intsvy)
mod3b <- intsvy::pirls.reg.pv( pvlabel="ASMMAT", x=c("migrant", "books" ),
             data=data.timss3a)

#--- regression with BIFIEsurvey
mod3c <- BIFIEsurvey::BIFIE.linreg( bdat4, dep="ASMMAT", pre=c("one","migrant","books"))
summary(mod3c)

#--- regression with lavaan.survey package
lavmodel <- "
    ASMMAT ~ migrant + books
    ASMMAT ~ 1
    ASMMAT ~~ ASMMAT
        "
# fit in lavaan
lavaan.fit <- lavaan::lavaan( lavmodel, data=data.timss3a, group="female", estimator="MLM")
summary(lavaan.fit)
mod3d <- lavaan.survey::lavaan.survey(lavaan.fit=lavaan.fit, survey.design=svydes4 )
summary(mod3d)

#############################################################################
# EXAMPLE 2: TIMSS dataset data.timss4 | Nested multiply imputed dataset
#############################################################################

data(data.timss4)
data(data.timssrep)

#**** create BIFIEdata object
bdat <- BIFIEsurvey::BIFIE.data( data.list=data.timss4, wgt=data.timss4[[1]][[1]]$TOTWGT,
               wgtrep=data.timssrep[, -1 ], NMI=TRUE, cdata=TRUE )
summary(bdat)

#**** Model 1: Linear regression for mathematics score
mod1 <- BIFIEsurvey::BIFIE.linreg( bdat, dep="ASMMAT", pre=c("one","books","migrant"))
summary(mod1)

#*** Model 2: Univariate statistics ?BIFIEsurvey::BIFIE.univar
mod2 <- BIFIEsurvey::BIFIE.univar( bdat, vars=c("ASMMAT","ASSSCI","books") )
summary(mod2)
}
}


\keyword{datasets}
