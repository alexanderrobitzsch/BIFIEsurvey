%% File Name: BIFIEsurvey-package.Rd
%% File Version: 0.9408

\name{BIFIEsurvey-package}
\alias{BIFIEsurvey-package}
\alias{BIFIEsurvey}
\docType{package}

\title{
Some Tools for Survey Statistics in Educational Assessment, Developed by BIFIE
}
\description{
    Contains tools for survey statistics (especially in educational
    assessment) for datasets with replication designs (jackknife, 
    bootstrap, replicate weights; see Kolenikov, 2010). 
    Descriptive statistics, linear and logistic regression, 
    path models for manifest variables with measurement error 
    correction and two-level hierarchical regressions for weighted 
    samples are included. Statistical inference can be conducted for 
    multiply imputed datasets and nested multiply imputed datasets
    and is in particularly suited for the analysis of plausible values
    (for details see George, Oberwimmer & Itzlinger-Bruneforth, 2016; 
    Bruneforth, Oberwimmer & Robitzsch, 2016; Robitzsch, Pham &
    Yanagida, 2016).    
    This package is developed by BIFIE (Federal Institute for 
    Educational Research, Innovation and Development of the Austrian 
    School System; Salzburg, Austria).
}

%% \code{Sys.time}: \tab syssys \cr

\details{

The \pkg{BIFIEsurvey} package include basic descriptive functions for large scale assessment data
to complement the more comprehensive \pkg{survey} package. The
functions in this package were written in \pkg{Rcpp}.

The features of \pkg{BIFIEsurvey} include for designs with replicate weights
(which includes Jackknife and Bootstrap as general approaches):

\itemize{

\item Descriptive statistics: means and standard deviations (\code{\link{BIFIE.univar}}),
frequencies (\code{\link{BIFIE.freq}}),
crosstabs (\code{\link{BIFIE.crosstab}})

\item Linear regression (\code{\link{BIFIE.linreg}})

\item Logistic regression (\code{\link{BIFIE.logistreg}})

\item Path models with measurement error correction for manifest
variables (\code{\link{BIFIE.pathmodel}})

\item Two-level regression for hierarchical data (\code{\link{BIFIE.twolevelreg}};
random slope model)

\item Statistical inference for derived parameters (\code{\link{BIFIE.derivedParameters}})

\item Wald tests (\code{\link{BIFIE.waldtest}}) of model parameters based on
replicated statistics

\item User-defined \R functions (\code{\link{BIFIE.by}})

}

}

%Statistical inference is conducted for multiply imputed datasets.
%Inference for nested multiple imputation (Rubin, 2003) is also included.


\author{
BIFIE [aut]

\url{http://www.bifie.at}

Bildungsforschung, Innovation & Entwicklung des oesterreichischen Schulwesens

Federal Institute for Educational Research, Innovation and
Development of the Austrian School System

Maintainer: Alexander Robitzsch <robitzsch@ipn.uni-kiel.de>

With contributions from: Konrad Oberwimmer, Alexander Robitzsch
}

%----
%For the following functions, an own function maintainer is defined
%which should be contacted instead of the package maintainer:
%
%Alexander Robitzsch <a.robitzsch@bifie.at> \cr
%\code{\link{BIFIE.logistreg}}
%----

\references{
Bruneforth, M., Oberwimmer, K., & Robitzsch, A. (2016). Reporting und Analysen. 
In S. Breit & C. Schreiner (Hrsg.). \emph{Large-Scale Assessment mit R: Methodische 
Grundlagen der oesterreichischen Bildungsstandardueberpruefung} (S. 333-362). Wien: facultas.

George, A. C., Oberwimmer, K., & Itzlinger-Bruneforth, U. (2016). Stichprobenziehung. 
In S. Breit & C. Schreiner (Hrsg.). \emph{Large-Scale Assessment mit R: Methodische 
Grundlagen der oesterreichischen Bildungsstandardueberpruefung} (S. 51-81). Wien: facultas.

Kolenikov, S. (2010). Resampling variance estimation for complex survey data.
\emph{Stata Journal, 10}(2), 165-199.

Robitzsch, A., Pham, G., & Yanagida, T. (2016). Fehlende Daten und Plausible Values. 
In S. Breit & C. Schreiner (Hrsg.). \emph{Large-Scale Assessment mit R: Methodische 
Grundlagen der oesterreichischen Bildungsstandardueberpruefung} (S. 259-293). Wien: facultas.
}

\keyword{package}
\seealso{
See also the \pkg{survey}, \pkg{intsvy}, \pkg{svyPVpack},
\pkg{lavaan.survey}, \pkg{EVER} and
the \pkg{eatRep} (\url{https://r-forge.r-project.org/R/?group_id=1326})
packages.
}

\examples{
##   |-----------------------------------------------------------------
##   | BIFIEsurvey 0.1-21 (2014-06-21)
##   | Maintainer: Alexander Robitzsch <a.robitzsch at bifie.at >
##   | http://www.bifie.at
##   |-----------------------------------------------------------------

##   .........................*,::;              :,:;       * .;*;. .,:
##   :::::::::::::::::::::::::.  ##+@              ##+#       .@####+ ;+# *
##   :::::::::::::::::::::::::.  ###@              ####       @@;   :*;##**
##   :::::::::::::::::::::::::.  ###@              ##+#     *##.      .,,
##   :::::::::::::::::::::::::.  ###@              ::,:    * ##
##   :::::::::::::::::::::::::.  ###@                 *      @@
##   :::::::::::::::::::::::::.  ###@                      * #@
##   :::::::::::::::::::::::::.  ###@                        #@                   *
##   :::::::::::::::::::::::::.  ##@#,@###@       @@## * @@@+##@@@@   #@     *.@#####
##   :::::::::::::::::::::::::.  ####*@#####@.**   ####   @@#@@@####   #@     ;##+**#+@*
##   :::::::::::::::::::::::::.  ##@@##,,#+##@     ####      @@        #@    .#@*  * *##
##   ::::,::::::::::::::::::::.  ##+@@    ####,    ####      #@        #@    ##,      #+*
##   ::::**.::::::::::::::::::.  ##+@, **,###@    ####      #@        #@    @@       ;@;
##   :::::*  .::::::::::::::::.  ##*       ####    ####      #@        #@   .##@@@@@@@#+#
##   :::::::*   *,,:::::::::::.     :      ###@    ####      #@        #@   ;##@@@@@@@@@@
##   :::::::::.      **....*    ,@#      ####    ####      #@        #@   *@**
##   ::::::::::::.*              ####      @###*   ####      #@        #@  *+#*
##   ::::::::::::::::,,,,,,,::.  ####      ####    ####      #@        #@  **#*
##   :::::::::::::::::::::::::.  ####      @###    ####      #@        #@ ,##
##   :::::::::::::::::::::::::.  ###@   *  @##+    ####      #@        #@   *#@
##   :::::::::::::::::::::::::*  @@##   ,#@##    ####      #@        #@    @#
##   :::::::::::::::::::::::::*  @+##,    @###*    ####      #@        #@   *+#@
##   ::::::::::::::::::::::::.    ###@:  *###@     ####      #@        #@     #@+
##   :::::::::::::::::::::::.   **;@#@#@####@.     ####      #@        #@     *@#@:*   *
##   ::::::::::::::::::::::.     ,@######@.      ####*     ##        @+      *#####@##
##   ::::::::::::::::::::.*        * .*##*. *        ***     *.        **         ;##+;.

}
