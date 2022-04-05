%% File Name: BIFIEsurvey-package.Rd
%% File Version: 1.021

\name{BIFIEsurvey-package}
\alias{BIFIEsurvey-package}
\alias{BIFIEsurvey}
\docType{package}

\title{
\packageTitle{BIFIEsurvey}
}

\description{
\packageDescription{BIFIEsurvey}
}

\author{
\packageAuthor{BIFIEsurvey}

Maintainer: \packageMaintainer{BIFIEsurvey}
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

Pfefferman, D., & Rao, C. R. (2009a). \emph{Handbook of statistics, Vol. 29A:
Sample surveys: Design, methods and applications}. Amsterdam: North Holland.

Pfefferman, D., & Rao, C. R. (2009b). \emph{Handbook of statistics, Vol. 29B:
Sample surveys: Inference and analysis}. Amsterdam: North Holland.

Robitzsch, A., Pham, G., & Yanagida, T. (2016). Fehlende Daten und Plausible Values. 
In S. Breit & C. Schreiner (Hrsg.). \emph{Large-Scale Assessment mit R: Methodische 
Grundlagen der oesterreichischen Bildungsstandardueberpruefung} (S. 259-293). Wien: facultas.

Shao, J. (1996). Invited discussion paper: Resampling methods in sample surveys. 
\emph{Statistics, 27}(3-4), 203-237.
}

\keyword{package}
\seealso{
See also the \pkg{survey}, \pkg{intsvy}, \pkg{EdSurvey},
\pkg{lavaan.survey}, \pkg{EVER} and
the \pkg{eatRep} packages.
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
