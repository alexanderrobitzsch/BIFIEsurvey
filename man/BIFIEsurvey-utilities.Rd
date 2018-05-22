%% File Name: BIFIEsurvey-utilities.Rd
%% File Version: 0.06

\name{BIFIEsurvey-utilities}
\alias{BIFIEsurvey-utilities}
\alias{bifiesurvey_rcpp_rubin_rules}
\alias{bifiesurvey_rcpp_replication_variance}
\alias{BIFIE_NMI_inference_parameters}


\title{Utility Functions in \pkg{BIFIEsurvey}}

\description{
Utility functions in \pkg{BIFIEsurvey}.
}

\usage{
## Rubin rules
bifiesurvey_rcpp_rubin_rules(estimates, variances)
## Replication variance
bifiesurvey_rcpp_replication_variance(pars, pars_repl, fay_factor)
## Inference for nested multiple imputation
BIFIE_NMI_inference_parameters( parsM, parsrepM, fayfac, RR, Nimp, Nimp_NMI,
      comp_cov=FALSE)
}

\arguments{
\item{estimates}{Vector}
\item{variances}{Vector}
\item{pars}{Matrix}
\item{pars_repl}{Matrix}
\item{fay_factor}{Vector}
\item{parsM}{Matrix}
\item{parsrepM}{Matrix}
\item{fayfac}{Vector}
\item{RR}{Numeric}
\item{Nimp}{Integer}
\item{Nimp_NMI}{Integer}
\item{comp_cov}{Logical}
}

%\details{
%xxx
%%}


% \examples{
% #############################################################################
% }
