%% File Name: BIFIEsurvey-utilities.Rd
%% File Version: 0.02
%% File Last Change: 2017-06-23 20:31:36

\name{BIFIEsurvey-utilities}
\alias{BIFIEsurvey-utilities}
\alias{bifie_rubin_rules}
\alias{bifie_replication_variance}
\alias{BIFIE_NMI_inference_parameters}


\title{Utility Functions in \pkg{BIFIEsurvey}}

\description{
Utility functions in \pkg{BIFIEsurvey}. 
}

\usage{
## Rubin rules
bifie_rubin_rules(estimates, variances)
## Replication variance
bifie_replication_variance(pars, pars_repl, fay_factor)
## Inference for nested multiple imputation 
BIFIE_NMI_inference_parameters( parsM, parsrepM, fayfac, RR, Nimp, Nimp_NMI, 
      comp_cov = FALSE) 
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