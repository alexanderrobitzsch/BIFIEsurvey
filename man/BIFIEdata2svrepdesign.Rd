%% File Name: BIFIEdata2svrepdesign.Rd
%% File Version: 0.02

\name{BIFIEdata2svrepdesign}
\alias{BIFIEdata2svrepdesign}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Conversion of a \code{BIFIEdata} Object into a \code{svyrep} Object in the
\pkg{survey} Package
}

\description{
Converts of a \code{BIFIEdata} object into a \code{svyrep} object in the
\pkg{survey} package.
}

\usage{
BIFIEdata2svrepdesign(bifieobj)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bifieobj}{
Object of class \code{BIFIEdata}
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
Object of class \code{svyrep.design} or \code{svyimputationList}
}


%%\references{
%% ~put references to the literature/web site here ~
%% }

%\author{
%%  ~~who you are~~
%}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See the \code{\link{BIFIE.data}} function for creating objects of class
\code{BIFIEdata} in \pkg{BIFIEsurvey}.

See the \code{\link[survey:svrepdesign]{survey::svrepdesign}} function in
the \pkg{survey} package.
}



\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: One dataset, TIMSS replication design
#############################################################################

data(data.timss3)
data(data.timssrep)

#--- create BIFIEdata object
bdat3 <- BIFIEsurvey::BIFIE.data.jack(data.timss3 , jktype="JK_TIMSS")
summary(bdat3)              

#--- create survey object directly in survey package
dat3a <- as.data.frame( cbind( data.timss3 , data.timssrep ) )
RR <- ncol(data.timssrep) - 1       # number of jackknife zones
svydes3a <- survey::svrepdesign(data=dat3a, weights=~TOTWGT ,type="JKn",
                 repweights='w_fstr[0-9]', scale=1,  rscales=rep(1,RR), mse=TRUE )
print(svydes3a)

#--- create survey object by converting the BIFIEdata object to survey
svydes3b <- BIFIEsurvey::BIFIEdata2svrepdesign(bdat3)

#*** compare results for the mean in Mathematics scores
mod1a <- BIFIEsurvey::BIFIE.univar( bdat3 , vars = "ASMMAT1")
mod1b <- survey::svymean( ~ ASMMAT1 , design=svydes3a )
mod1c <- survey::svymean( ~ ASMMAT1 , design=svydes3b )

#- coefficients
coef(mod1a); coef(mod1b); coef(mod1c)
#- standard errors
survey::SE(mod1a); survey::SE(mod1b); survey::SE(mod1c)

#############################################################################
# EXAMPLE 2: Multiply imputed datasets, TIMSS replication design
#############################################################################

data(data.timss2)
data(data.timssrep)

#--- create BIFIEdata object
bdat4 <- BIFIEsurvey::BIFIE.data( data=data.timss2, wgt = "TOTWGT", 
              wgtrep= data.timssrep[,-1], fayfac = 1)
print(bdat4)

#--- create object with imputed datasets in survey
datL <- mitools::imputationList( data.timss2 )    
RR <- ncol(data.timssrep) - 1
weights <- data.timss2[[1]]$TOTWGT
repweights <-  data.timssrep[,-1]
svydes4a <- survey::svrepdesign(data=datL, weights= weights, type="other", 
               repweights= repweights , scale=1,  rscales=rep(1,RR), mse=TRUE)
print(svydes4a)

#--- create BIFIEdata object with conversion function
svydes4b <- BIFIEsurvey::BIFIEdata2svrepdesign(bdat4)

#*** compare results for a mean
mod1a <- BIFIEsurvey::BIFIE.univar(bdat4 , vars = "ASMMAT")
mod1b <- mitools::MIcombine( with(svydes4a, survey::svymean( ~ ASMMAT, design=svydes4a )))
mod1c <- mitools::MIcombine( with(svydes4b, survey::svymean( ~ ASMMAT, design=svydes4b )))
# results
coef(mod1a); coef(mod1b); coef(mod1c)
survey::SE(mod1a); survey::SE(mod1b); survey::SE(mod1c)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
