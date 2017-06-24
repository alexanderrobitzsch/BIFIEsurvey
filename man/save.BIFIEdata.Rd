\name{save.BIFIEdata}
\alias{save.BIFIEdata}
\alias{load.BIFIEdata}
\alias{write.BIFIEdata}
\alias{load.BIFIEdata.files}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Saving, Writing and Loading \code{BIFIEdata} Objects
}

\description{
These functions save (\code{save.BIFIEdata}) , write (\code{write.BIFIEdata})
or load (\code{load.BIFIEdata}) objects of class \code{BIFIEdata}.

The function \code{load.BIFIEdata.files} allows the creation
of \code{BIFIEdata} objects by loading separate files of imputed datasets,
replicate weights and a possible indicator dataset.
}

\usage{
save.BIFIEdata(BIFIEdata, name.BIFIEdata , cdata=TRUE , varnames=NULL)

write.BIFIEdata( BIFIEdata , name.BIFIEdata , dir = getwd() , varnames = NULL ,
    impdata.index = NULL , type = "Rdata" , ... )

load.BIFIEdata(filename, dir=getwd() )

load.BIFIEdata.files( files.imp , wgt , file.wgtrep , file.ind = NULL , 
    type="Rdata" ,varnames=NULL , cdata =TRUE , dir=getwd() , ... )
}

%% label the arguments of load.BIFIEdata.files
%- maybe also 'usage' for other objects documented here.


\arguments{
  \item{BIFIEdata}{
Object of class \code{BIFIEdata}
}
\item{name.BIFIEdata}{
	Name of \code{BIFIEdata} set to be saved
   }
\item{cdata}{An optional logical indicating whether the dataset should be
saved in a 'compact way'}
\item{varnames}{Vector of variable names which should be saved. The default
	is to use all variables.	
}
\item{dir}{Directory in which data files should be saved. The default
is the working directory.}
\item{impdata.index}{Vector of indices for selecting imputed datasets}
\item{type}{Type of saved data. Options are \code{Rdata} (function 
	\code{\link[base:save]{base::save}}, 
	\code{csv} (function \code{\link[utils:write.csv]{utils::write.csv}}),
	\code{csv2} (function \code{\link[utils:write.csv2]{utils::write.csv2}}),
	\code{table} (function \code{\link[utils:write.table]{utils::write.table}}),
	\code{sav} (function \code{\link[foreign:read.spss]{foreign::read.spss}}
		for reading sav files and function
	\code{\link[sjmisc:write_spss]{sjmisc::write_spss}}.
		}
\item{\dots}{Additional arguments to be passed to 
	\code{\link[base:save]{base::save}},
	\code{\link[utils:write.csv]{utils::write.csv}},
	\code{\link[utils:write.csv2]{utils::write.csv2}},
	\code{\link[utils:write.table]{utils::write.table}}, 
	\code{\link[foreign:read.spss]{foreign::read.spss}},
	\code{\link[sjmisc:write_spss]{sjmisc::write_spss}}
}
\item{filename}{
	File name of \code{BIFIEdata} object
   }
\item{files.imp}{Vector of file names of imputed datasets}
\item{wgt}{Variable name of case weight}
\item{file.wgtrep}{File name for dataset with replicate weights}
\item{file.ind}{Optional. File name for dataset with response data indicators}
}

% load.BIFIEdata.files( files.imp , wgt , file.wgtrep , file.ind = NULL , 
%    type="Rdata" ,varnames=NULL , cdata =TRUE , dir=getwd() , ... )


%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
Saved \R object and a summary in working directory or a loaded \R object. 
}

%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%For saving files in sav format, an installed version of the freeware
%software PSPP is necessary, see \code{\link[miceadds:write.pspp]{write.pspp}}.
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For creating objects of class \code{BIFIEdata} see \code{\link{BIFIE.data}}.

\code{\link[base:save]{base::save}}, \code{\link[base:load]{base::load}}
}



\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Saving and loading BIFIE data objects
#############################################################################
data(data.timss1)
data(data.timssrep)

bifieobj <- BIFIEsurvey::BIFIE.data( data.list=data.timss1, wgt=data.timss1[[1]]$TOTWGT,
                 wgtrep=data.timssrep[, -1 ] )
summary(bifieobj)       

# save bifieobj in a compact way 
save.BIFIEdata( BIFIEdata=bifieobj, name.BIFIEdata="timss1_cdata" )
# save bifieobj in a non-compact way
save.BIFIEdata( BIFIEdata=bifieobj, name.BIFIEdata="timss1_data", cdata=FALSE )

# load this object with object name "bdat2"
bdat2 <- BIFIEsurvey::load.BIFIEdata( filename = "timss1_data.Rdata" )
summary(bdat2)

# save bifieobj with selected variables
save.BIFIEdata( bifieobj, name.BIFIEdata="timss1_selectvars_cdata" , 
                     varnames= bifieobj$varnames[ c(1:7,13,12,9) ] )
# the same object, but use the non-compact way of saving
save.BIFIEdata( bifieobj, name.BIFIEdata="timss1_selectvars_data", cdata=FALSE, 
                     varnames= bifieobj$varnames[ c(1:7,13,12,9) ] )

# load object timss1_cdata (in compact data format)
bdat3 <- BIFIEsurvey::load.BIFIEdata( filename = "timss1_cdata.Rdata" )
summary(bdat3)
# save selected variables of object bdat3
save.BIFIEdata( bdat3, name.BIFIEdata="timss1_selectvars2_cdata" ,  
                     varnames= bifieobj$varnames[ c(1:4,12,8) ] )      
                     
#############################################################################
# EXAMPLE 2: Writing BIFIEdata objects 
#############################################################################                     

data(data.timss2)
data(data.timssrep)

# create compactBIFIEdata
bifieobj <- BIFIEsurvey::BIFIE.data( data.list=data.timss2, wgt=data.timss2[[1]]$TOTWGT,
                wgtrep=data.timssrep[, -1 ] , cdata=TRUE)
summary(bifieobj)           

# save imputed datasets in format csv2
write.BIFIEdata( bifieobj, name.BIFIEdata="timss2_save1", type="csv2", row.names=FALSE)

# save imputed datasets of BIFIEdata object in format table without column names 
# and code missings as "."
write.BIFIEdata( bifieobj , name.BIFIEdata = "timss2_save2" , type="table" , 
                  col.names=FALSE , row.names=FALSE , na="." )

# save imputed datasets of  BIFIEdata object in format csv and select some variables 
# and only the first three datasets
varnames <- c("IDSTUD","TOTWGT","female","books","lang","ASMMAT")
write.BIFIEdata( bifieobj , name.BIFIEdata = "timss2_save3" , type="csv" ,
                   impdata.index = 1:3  , varnames=varnames)

# save imputed datasets of BIFIEdata object in format Rdata, the R binary format            
write.BIFIEdata( bifieobj , name.BIFIEdata = "timss2_save4" , type="Rdata"  )

# save imputed datasets in sav (SPSS) format
write.BIFIEdata( bifieobj , name.BIFIEdata = "timss2_save5" , type="sav" )
             
#############################################################################
# EXAMPLE 3: Loading BIFIEdata objects saved in separate files 
#                   (no indicator dataset)
#############################################################################              

# We assume that Example 2 is applied and we build on the saved files
# from this example.

#***--- read Rdata format
# extract files with imputed datasets and replicate weights
files.imp <- miceadds::grep.vec( c("timss2_save4__IMP" , ".Rdata" ) , 
			list.files(getwd())  )$x
file.wgtrep <- miceadds::grep.vec( c("timss2_save4__WGTREP" , ".Rdata" ) , 
			list.files(getwd())  )$x
# select some variables in varnames
varnames <- scan( nlines=1, what="character")
   IDSTUD   TOTWGT books lang migrant likesc  ASMMAT

# load files and create BIFIEdata object
bifieobj1 <- BIFIEsurvey::load.BIFIEdata.files( files.imp, wgt="TOTWGT", file.wgtrep,
                        type="Rdata" , varnames=varnames )
summary(bifieobj1)

#***--- read csv2 format
files.imp <- miceadds::grep.vec( c("timss2_save1__IMP" , ".csv" ) ,
                        list.files(getwd()) )$x
file.wgtrep <- miceadds::grep.vec( c("timss2_save1__WGTREP" , ".csv" ) , 
                        list.files(getwd()) )$x
bifieobj2 <- BIFIEsurvey::load.BIFIEdata.files( files.imp, wgt="TOTWGT", 
                    file.wgtrep , type="csv2" )
summary(bifieobj2)

#***--- read sav format
files.imp <- miceadds::grep.vec( c("timss2_save5__IMP" , ".sav" ) , 
                        list.files(getwd()) )$x
file.wgtrep <- miceadds::grep.vec( c("timss2_save5__WGTREP" , ".sav" ) , 
                        list.files(getwd()) )$x
bifieobj3 <- BIFIEsurvey::load.BIFIEdata.files( files.imp , wgt="TOTWGT", 
                file.wgtrep, type="sav", to.data.frame=TRUE , use.value.labels=FALSE)
summary(bifieobj3)

#############################################################################
# EXAMPLE 4: Loading BIFIEdata objects saved in separate files 
#                   (with an indicator dataset)
############################################################################# 

data(data.timss1)
data(data.timss1.ind)
data(data.timssrep)

# create BIFIEdata object at first
bifieobj <- BIFIEsurvey::BIFIE.data( data.list=data.timss1 , wgt= "TOTWGT" ,
            wgtrep=data.timssrep[, -1 ] )
summary(bifieobj)       

#--- save datasets for the purpose of the following example
write.BIFIEdata( BIFIEdata=bifieobj, name.BIFIEdata="timss1_ex" , type="Rdata" )
# save indicator dataset
save( data.timss1.ind , file="timss1_ex__IND.Rdata" )

# grep file names
files.imp <- miceadds::grep.vec( c("timss1_ex__IMP" , ".Rdata" ) , 
                      list.files(getwd()) )$x
file.wgtrep <- miceadds::grep.vec( c("timss1_ex__WGTREP" , ".Rdata" ) , 
                      list.files(getwd()) )$x
file.ind <- miceadds::grep.vec( c("timss1_ex__IND" , ".Rdata" ) , 
                      list.files(getwd()) )$x
# define variables for selection
varnames <- c("IDSTUD","TOTWGT","female","books","lang","ASMMAT")
# read files using indicator dataset
bifieobj2 <- BIFIEsurvey::load.BIFIEdata.files( files.imp , wgt="TOTWGT", 
                  file.wgtrep=file.wgtrep, file.ind=file.ind, type="Rdata", 
                  varnames=varnames)
summary(bifieobj2)

# read files without indicator dataset
bifieobj3 <- BIFIEsurvey::load.BIFIEdata.files( files.imp , wgt="TOTWGT", 
                file.wgtrep=file.wgtrep, type="Rdata"  , varnames = varnames)
summary(bifieobj3)

# compare some descriptive statistics
res2 <- BIFIEsurvey::BIFIE.univar( bifieobj2 , vars = c("books" , "ASMMAT" , "lang") )
res3 <- BIFIEsurvey::BIFIE.univar( bifieobj3 , vars = c("books" , "ASMMAT" , "lang") )
summary(res2)
summary(res3)
}                     
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{BIFIEdata}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
