//// File Name: bifiesurvey_rcpp_wald_test.h
//// File Version: 0.175

#ifndef _BIFIESURVEY_BIFIESURVEY_RCPP_WALD_TEST_H
#define _BIFIESURVEY_BIFIESURVEY_RCPP_WALD_TEST_H
 
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
#include "bifiesurvey_rcpp_helper.h"
using namespace Rcpp;
using namespace arma;

Rcpp::List bifiesurvey_rcpp_wald_test_vcov( int VV, Rcpp::NumericVector Ccols,
    Rcpp::NumericMatrix parsM, Rcpp::NumericMatrix parsrepM, int ii, int RR,
    Rcpp::NumericVector fayfac, arma::mat ACdes, arma::colvec Ardes );

Rcpp::List bifiesurvey_rcpp_wald_test( Rcpp::NumericMatrix parsM, Rcpp::NumericMatrix parsrepM,
    Rcpp::NumericMatrix Cdes, Rcpp::NumericVector rdes, Rcpp::NumericVector Ccols,
    Rcpp::NumericVector fayfac );

#endif // _BIFIESURVEY_BIFIESURVEY_RCPP_WALD_TEST_H
