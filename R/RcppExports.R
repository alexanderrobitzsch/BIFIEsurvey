## File Name: RcppExports.R
## File Version: 3.007001
# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

bifiesurvey_rcpp_jackknife_timss <- function(wgt, jkzone, jkrep, RR, jkfac, prbar) {
    .Call('_BIFIEsurvey_bifiesurvey_rcpp_jackknife_timss', PACKAGE='BIFIEsurvey', wgt, jkzone, jkrep, RR, jkfac, prbar)
}

bifiesurvey_rcpp_bootstrap <- function(cumwgt, rand_wgt) {
    .Call('_BIFIEsurvey_bifiesurvey_rcpp_bootstrap', PACKAGE='BIFIEsurvey', cumwgt, rand_wgt)
}

bifiesurvey_rcpp_bifiedata2bifiecdata <- function(datalistM, Nimp) {
    .Call('_BIFIEsurvey_bifiesurvey_rcpp_bifiedata2bifiecdata', PACKAGE='BIFIEsurvey', datalistM, Nimp)
}

bifiesurvey_rcpp_bifiecdata2bifiedata <- function(datalistM_ind, datalistM_imputed, Nimp, dat1, datalistM_impindex) {
    .Call('_BIFIEsurvey_bifiesurvey_rcpp_bifiecdata2bifiedata', PACKAGE='BIFIEsurvey', datalistM_ind, datalistM_imputed, Nimp, dat1, datalistM_impindex)
}

bifiesurvey_rcpp_bifiedata_stepwise <- function(dat1, dat_ind, Nmiss) {
    .Call('_BIFIEsurvey_bifiesurvey_rcpp_bifiedata_stepwise', PACKAGE='BIFIEsurvey', dat1, dat_ind, Nmiss)
}

bifiesurvey_rcpp_linreg <- function(datalist, wgt1, wgtrep, dep_index, pre_index, fayfac, NI, group_index1, group_values) {
    .Call('_BIFIEsurvey_bifiesurvey_rcpp_linreg', PACKAGE='BIFIEsurvey', datalist, wgt1, wgtrep, dep_index, pre_index, fayfac, NI, group_index1, group_values)
}

bifiesurvey_rcpp_logistreg <- function(datalist, wgt1, wgtrep, dep_index, pre_index, fayfac, NI, group_index1, group_values, eps, maxiter) {
    .Call('_BIFIEsurvey_bifiesurvey_rcpp_logistreg', PACKAGE='BIFIEsurvey', datalist, wgt1, wgtrep, dep_index, pre_index, fayfac, NI, group_index1, group_values, eps, maxiter)
}

univar_multiple_V2group <- function(datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values) {
    .Call('_BIFIEsurvey_univar_multiple_V2group', PACKAGE='BIFIEsurvey', datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values)
}

bifie_freq <- function(datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values, vars_values, vars_values_numb) {
    .Call('_BIFIEsurvey_bifie_freq', PACKAGE='BIFIEsurvey', datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values, vars_values, vars_values_numb)
}

bifie_correl <- function(datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values) {
    .Call('_BIFIEsurvey_bifie_correl', PACKAGE='BIFIEsurvey', datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values)
}

bifie_comp_vcov_within <- function(parsM, parsrepM, fayfac, RR, Nimp) {
    .Call('_BIFIEsurvey_bifie_comp_vcov_within', PACKAGE='BIFIEsurvey', parsM, parsrepM, fayfac, RR, Nimp)
}

bifie_comp_vcov <- function(parsM, parsrepM, Cdes, rdes, Ccols, fayfac) {
    .Call('_BIFIEsurvey_bifie_comp_vcov', PACKAGE='BIFIEsurvey', parsM, parsrepM, Cdes, rdes, Ccols, fayfac)
}

bifie_test_univar <- function(mean1M, sd1M, sumweightM, GG, group_values, mean1repM, sd1repM, sumweightrepM, fayfac) {
    .Call('_BIFIEsurvey_bifie_test_univar', PACKAGE='BIFIEsurvey', mean1M, sd1M, sumweightM, GG, group_values, mean1repM, sd1repM, sumweightrepM, fayfac)
}

bifie_crosstab <- function(datalist, wgt1, wgtrep, vars_values1, vars_index1, vars_values2, vars_index2, fayfac, NI, group_index1, group_values) {
    .Call('_BIFIEsurvey_bifie_crosstab', PACKAGE='BIFIEsurvey', datalist, wgt1, wgtrep, vars_values1, vars_index1, vars_values2, vars_index2, fayfac, NI, group_index1, group_values)
}

bifie_by <- function(datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values, userfct) {
    .Call('_BIFIEsurvey_bifie_by', PACKAGE='BIFIEsurvey', datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values, userfct)
}

bifie_hist <- function(datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values, breaks) {
    .Call('_BIFIEsurvey_bifie_hist', PACKAGE='BIFIEsurvey', datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values, breaks)
}

bifie_ecdf <- function(datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values, breaks, quanttype, maxval) {
    .Call('_BIFIEsurvey_bifie_ecdf', PACKAGE='BIFIEsurvey', datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values, breaks, quanttype, maxval)
}

bifie_fasttable <- function(datavec) {
    .Call('_BIFIEsurvey_bifie_fasttable', PACKAGE='BIFIEsurvey', datavec)
}

bifie_table1_character <- function(datavec) {
    .Call('_BIFIEsurvey_bifie_table1_character', PACKAGE='BIFIEsurvey', datavec)
}

bifie_mla2 <- function(X_list, Z_list, y_list, wgttot, wgtlev2, wgtlev1, globconv, maxiter, group, group_values, cluster, wgtrep, Nimp, fayfac, recov_constraint, is_rcov_constraint) {
    .Call('_BIFIEsurvey_bifie_mla2', PACKAGE='BIFIEsurvey', X_list, Z_list, y_list, wgttot, wgtlev2, wgtlev1, globconv, maxiter, group, group_values, cluster, wgtrep, Nimp, fayfac, recov_constraint, is_rcov_constraint)
}

bifiesurvey_rcpp_replication_variance <- function(pars, pars_repl, fay_factor) {
    .Call('_BIFIEsurvey_bifiesurvey_rcpp_replication_variance', PACKAGE='BIFIEsurvey', pars, pars_repl, fay_factor)
}

bifiesurvey_rcpp_rubin_rules <- function(estimates, variances) {
    .Call('_BIFIEsurvey_bifiesurvey_rcpp_rubin_rules', PACKAGE='BIFIEsurvey', estimates, variances)
}

bifiesurvey_rcpp_pathmodel <- function(datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values, L, L_row_index, NL, E, R, R_row_index, coeff_index, NP0, unreliability) {
    .Call('_BIFIEsurvey_bifiesurvey_rcpp_pathmodel', PACKAGE='BIFIEsurvey', datalist, wgt1, wgtrep, vars_index, fayfac, NI, group_index1, group_values, L, L_row_index, NL, E, R, R_row_index, coeff_index, NP0, unreliability)
}

bifiesurvey_rcpp_wald_test <- function(parsM, parsrepM, Cdes, rdes, Ccols, fayfac) {
    .Call('_BIFIEsurvey_bifiesurvey_rcpp_wald_test', PACKAGE='BIFIEsurvey', parsM, parsrepM, Cdes, rdes, Ccols, fayfac)
}

