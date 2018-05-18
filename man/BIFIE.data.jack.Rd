%% File Name: BIFIE.data.jack.Rd
%% File Version: 0.26

\name{BIFIE.data.jack}
\alias{BIFIE.data.jack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create \code{BIFIE.data} Object with Jackknife Zones
}

\description{
Creates a \code{BIFIE.data} object for designs with jackknife zones,
especially for TIMSS/PIRLS and PISA studies.
}

\usage{
BIFIE.data.jack(data, wgt = NULL, jktype = "JK_TIMSS", pv_vars = NULL, 
     jkzone = NULL, jkrep = NULL, jkfac = NULL, fayfac = NULL , 
     wgtrep = "W_FSTR" , pvpre = paste0("PV",1:5) , ngr=100 ,
     seed = .Random.seed , cdata=FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data frame: Can be a single or a list of multiply-imputed datasets
}
  \item{wgt}{
A string indicating the label of case weight.
In case of \code{jktype="JK_TIMSS"} the weight is specified as
 \code{wgt="TOTWGT"} as the default.
}
  \item{pv_vars}{
An optional vector of plausible values which define multiply-imputed datasets.
}
  \item{jktype}{
Type of jackknife procedure for creating the \code{BIFIE.data} object.
\code{jktype="JK_TIMSS"} refers to TIMSS/PIRLS datasets up to 2011 data,
\code{jktype="JK_TIMSS2"} refers to TIMSS/PIRLS datasets starting from 2015 data. The type
\code{"JK_GROUP"} creates jackknife weights based on a user defined
grouping, the type \code{"JK_RANDOM"} creates random groups. The number
of random groups can be defined in \code{ngr}. The argument \code{type="RW_PISA"}
converts PISA datasets into objects of class \code{BIFIEdata}.
}
  \item{jkzone}{
Jackknife zones. If \code{jktype="JK_TIMSS"}, then \code{jkzone="JKZONE"}.
}
  \item{jkrep}{
Jackknife replicate factors. If \code{jktype="JK_TIMSS"}, then \code{jkrep="JKREP"}.
}
  \item{jkfac}{
Factor for multiplying jackknife replicate weights. 
If \code{jktype="JK_TIMSS"}, then \code{jkfac=2}.
}
  \item{fayfac}{
Fay factor for statistical inference. The default is set to \code{NULL}.
}
\item{wgtrep}{Variables in the dataset which refer to the replicate
    weights. In case of \code{cdata=TRUE}, the replicate weights
    are deleted from \code{datalistM}.
        }
\item{pvpre}{Only applicable for \code{jktype="RW_PISA"}. The vector contains
    the prefixes of the variables containing plausible values.    
}
\item{ngr}{Number of randomly created groups in \code{"JK_RANDOM"}.}
\item{seed}{The simulation seed if \code{"JK_RANDOM"} is chosen.
    If \code{seed=NULL}, then the grouping is done according the
order in the dataset.}
\item{cdata}{An optional logical indicating whether the \code{BIFIEdata}
object should be compactly saved. The default is \code{FALSE}.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
Object of class \code{BIFIEdata}
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
\code{\link{BIFIE.data}}, \code{\link{BIFIE.data.boot}}
}



\examples{
#############################################################################
# EXAMPLE 1: Convert TIMSS dataset to BIFIE.data object
#############################################################################

data(data.timss3)

# define plausible values
pv_vars <- c("ASMMAT" , "ASSSCI" )
# create BIFIE.data objects -> 5 imputed datasets
bdat1 <- BIFIEsurvey::BIFIE.data.jack( data=data.timss3  ,  pv_vars = pv_vars ,
             jktype="JK_TIMSS"  )
summary(bdat1)

# create BIFIE.data objects -> all PVs are included in one dataset
bdat2 <- BIFIEsurvey::BIFIE.data.jack( data=data.timss3  ,  jktype="JK_TIMSS"  )
summary(bdat2)

#############################################################################
# EXAMPLE 2: Creation of Jackknife zones and replicate weights for data.test1
#############################################################################

data(data.test1)

# create jackknife zones based on random group creation
bdat1 <- BIFIEsurvey::BIFIE.data.jack( data = data.test1 ,  jktype="JK_RANDOM", 
                    ngr = 50 )
summary(bdat1)
stat1 <- BIFIEsurvey::BIFIE.univar( bdat1 , vars ="math" ,  group="stratum" )
summary(stat1)

# random creation of groups and inclusion of weights
bdat2 <- BIFIEsurvey::BIFIE.data.jack( data = data.test1 ,  jktype="JK_RANDOM", 
                ngr = 75, seed=987, wgt="wgtstud")
summary(bdat2)
stat2 <- BIFIEsurvey::BIFIE.univar( bdat2 , vars ="math" ,  group="stratum" )
summary(stat2)

# using idclass as jackknife zones
bdat3 <- BIFIEsurvey::BIFIE.data.jack( data = data.test1 ,  jktype="JK_GROUP", 
                jkzone="idclass", wgt="wgtstud")
summary(bdat3)
stat3 <- BIFIEsurvey::BIFIE.univar( bdat3 , vars ="math" ,  group="stratum" )
summary(stat3)

# create BIFIEdata object with a list of imputed datasets
dataList <- list( data.test1 , data.test1 , data.test1 )
bdat4 <- BIFIEsurvey::BIFIE.data.jack( data = dataList ,  jktype="JK_GROUP" , 
                jkzone="idclass", wgt="wgtstud")
summary(bdat4)                

\dontrun{
#############################################################################
# EXAMPLE 3: Converting a PISA dataset into a BIFIEdata object
#############################################################################

data(data.pisaNLD)

# BIFIEdata with cdata=FALSE
bifieobj <- BIFIEsurvey::BIFIE.data.jack( data.pisaNLD, jktype="RW_PISA", cdata=FALSE)
summary(bifieobj)
# BIFIEdata with cdata=TRUE
bifieobj1 <- BIFIEsurvey::BIFIE.data.jack( data.pisaNLD, jktype="RW_PISA", cdata=TRUE)
summary(bifieobj1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{BIFIEdata}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
