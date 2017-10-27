%% File Name: data.pisaNLD.Rd
%% File Version: 0.13

\name{data.pisaNLD}
\alias{data.pisaNLD}

\docType{data}
\title{
Some PISA Datatsets
}

\description{
Some PISA datasets.
}

\usage{
data(data.pisaNLD)
}


\format{
The dataset \code{data.pisaNLD} is a data frame with 3992 observations on 405 variables
which is a part of the Dutch PISA 2006 data.
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Downloaded from \url{http://www.jstatsoft.org/v20/i05/}
}


%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
\dontrun{	
library(mitools)
library(survey)
library(intsvy)

#############################################################################
# EXAMPLE 1: Dutch PISA 2006 dataset
#############################################################################  

data(data.pisaNLD)
data <- data.pisaNLD

#--- Create object of class BIFIEdata

# list variables with plausible values: These must be named
# as pv1math, pv2math, ... , pv5math, ...
pv_vars <- toupper( c("math" , "math1" , "math2" , "math3" , "math4" ,
             "read" , "scie" , "prob") )
# create 5 datasets including different sets of plausible values
dfr <- NULL
VV <- length(pv_vars)
Nimp <- 5           # number of plausible values
for (vv in 1:VV){
      vv1 <- pv_vars[vv]
      ind.vv1 <- which( colnames(data) \%in\% paste0("PV" , 1:Nimp , vv1) )      
      dfr2 <- data.frame( "variable" = paste0("PV" , vv1) , "var_index" = vv , 
          "data_index" = ind.vv1 , "impdata_index"=1:Nimp )           
      dfr <- rbind( dfr , dfr2 )
}
       
sel_ind <- setdiff( 1:( ncol(data) ) , dfr$data_index )
data0 <- data[ , sel_ind ]
V0 <- ncol(data0)
newvars <- seq( V0+1 , V0+VV )
datalist <- as.list( 1:Nimp )
for (ii in 1:Nimp ){
    dat1 <- data.frame( data0 , data[ , dfr[ dfr$impdata_index == ii , "data_index" ]])
    colnames(dat1)[ newvars ] <- paste0("PV",pv_vars)
    datalist[[ii]] <- dat1 
}

# dataset with replicate weights
datarep <- data[ , grep( "W_FSTR" , colnames(data) ) ]
RR <- ncol(datarep)     # number of replicate weights

# create BIFIE object
bifieobj <- BIFIEsurvey::BIFIE.data( datalist , wgt = data[, "W_FSTUWT"], 
                 wgtrep = datarep, fayfac = 1 / RR / ( 1 - .5 )^2 )
# For PISA: RR=80 and therefore fayfac = 1/20 = .05
summary(bifieobj)

#--- Create BIFIEdata object immediately using BIFIE.data.jack function
bifieobj1 <- BIFIEsurvey::BIFIE.data.jack( data.pisaNLD, jktype = "RW_PISA", cdata=TRUE)
summary(bifieobj1)

#--- Create object in survey package
datL <- mitools::imputationList(list( datalist[[1]],datalist[[2]],
                  datalist[[3]],datalist[[4]],datalist[[5]]) )
pisades <- survey::svrepdesign(ids = ~ 1, weights = ~W_FSTUWT, data = datL ,
                    repweights = "W_FSTR[0-9]+", type = "Fay", rho = 0.5, mse=TRUE)
print(pisades)

#++++++++++++++ some comparisons with other packages +++++++++++++++++++++++++++++++

#**** Model 1: Means for mathematics and reading
# BIFIEsurvey package
mod1a <- BIFIEsurvey::BIFIE.univar( bifieobj , vars = c("PVMATH" , "PVREAD") )
summary(mod1a)

# intsvy package
mod1b <- intsvy::pisa.mean.pv(pvlabel = "MATH", data = data.pisaNLD )
mod1b

# survey package
mod1c <- with( pisades , survey::svymean(PVMATH~1, design=pisades) )
res1c <- mitools::MIcombine(mod1c)
summary(res1c)

#**** Model 2: Linear regression
# BIFIEsurvey package
mod2a <- BIFIEsurvey::BIFIE.linreg( bifieobj, dep="PVMATH", 
              pre=c("one","ANXMAT","HISEI"))
summary(mod2a)

# intsvy package
mod2b <- intsvy::pisa.reg.pv(pvlabel="MATH", x= c("ANXMAT","HISEI"), data=data.pisaNLD)
mod2b

# survey package
mod2c <- with( pisades , survey::svyglm(PVMATH~ANXMAT+HISEI, design=pisades) )
res2c <- mitools::MIcombine(mod2c)
summary(res2c)
}
}


\keyword{datasets}
