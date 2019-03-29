//// File Name: bifiesurvey_rcpp_helper.h
//// File Version: 7.502

#ifndef _BIFIESURVEY_BIFIESURVEY_RCPP_HELPER_H
#define _BIFIESURVEY_BIFIESURVEY_RCPP_HELPER_H
 
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double bifiesurvey_rcpp_squeeze( double x, double min_val, double max_val );

double bifiesurvey_rcpp_arma_trace( arma::mat x );

double bifiesurvey_rcpp_extract_fayfac( Rcpp::NumericVector fayfac, int rr);

Rcpp::NumericVector bifie_sign( double x );

Rcpp::NumericVector matr2vec( Rcpp::NumericMatrix matr1);

Rcpp::List matrix_entry( Rcpp::NumericMatrix parsM, Rcpp::NumericMatrix pars_full1,
        int vv );

Rcpp::List univar_helper_multiple_V2group( Rcpp::NumericMatrix dat1,
     Rcpp::NumericMatrix wgt1,  Rcpp::NumericVector vars_index,
     Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values );

Rcpp::NumericVector varjack_helper( Rcpp::NumericVector pars,
    Rcpp::NumericMatrix pars_jack, Rcpp::NumericVector fayfac );

Rcpp::List varjack_bias_helper( Rcpp::NumericVector pars,
    Rcpp::NumericMatrix pars_jack, Rcpp::NumericVector fayfac );

Rcpp::List rubin_rules_univ( Rcpp::NumericMatrix parsM, Rcpp::NumericMatrix pars_varM);

Rcpp::List bifiehelper_freq( Rcpp::NumericMatrix dat1, Rcpp::NumericMatrix wgt,
    Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values,
    Rcpp::NumericVector vars_values_numb, Rcpp::NumericMatrix vars_values,
    Rcpp::NumericVector vars_index, Rcpp::NumericVector vars_values_numb_cumsum );

Rcpp::List bifiehelpers_correl( Rcpp::NumericMatrix dat1, Rcpp::NumericVector ind_cases,
    Rcpp::NumericVector group_values, Rcpp::NumericVector group_index1,
    Rcpp::NumericMatrix wgt, Rcpp::NumericVector vars_index,
    Rcpp::NumericMatrix itempair_index);

Rcpp::List bifiehelpers_etasquared( Rcpp::NumericMatrix mean1M,
    Rcpp::NumericMatrix sd1M, Rcpp::NumericMatrix sumweightM, int GG );

Rcpp::List bifiehelpers_crosstab( Rcpp::NumericMatrix dat1, Rcpp::NumericMatrix wgt,
    Rcpp::NumericVector group_values, Rcpp::NumericVector group_index1,
    Rcpp::NumericVector vars_values1, Rcpp::NumericVector vars_index1,
    Rcpp::NumericVector vars_values2, Rcpp::NumericVector vars_index2,
    Rcpp::NumericMatrix design_pars );

Rcpp::NumericVector bifie_helper_ecdf( Rcpp::NumericMatrix dat1,
    Rcpp::NumericVector wgt, Rcpp::NumericVector breaks,
    Rcpp::NumericVector group_values, Rcpp::NumericVector group_index1,
    Rcpp::NumericVector vars_index, int ii,
    Rcpp::NumericMatrix ncasesM, Rcpp::NumericMatrix sumwgtM,
    int maxval, int quanttype );

Rcpp::List rcppmat2armamat( Rcpp::NumericMatrix matr );

Rcpp::List mla2_se_fixed( arma::mat theta, arma::mat Tmat, arma::mat sig2,
        arma::mat Var_theta_rj, arma::mat AfjArj,
        Rcpp::NumericVector wgtlev2, arma::mat Xt2 );

Rcpp::List mla2_decomp( Rcpp::NumericMatrix V,
    Rcpp::NumericMatrix idcluster_table, Rcpp::NumericVector wgttot_ );

Rcpp::NumericVector maxabsval_arma( arma::mat Tmat, arma::mat Tmat0 );

Rcpp::NumericVector mla2_checkconv( arma::mat theta, arma::mat theta0,
    arma::mat Tmat, arma::mat Tmat0, arma::mat sig2,
    arma::mat sig20 );

Rcpp::NumericVector mla2_vardec( arma::mat theta, arma::mat Tmat, arma::mat sig2,
    Rcpp::NumericMatrix Sigma_B_yX,Rcpp::NumericMatrix Sigma_W_yX,
    Rcpp::NumericMatrix Sigma_B_yZ,Rcpp::NumericMatrix Sigma_W_yZ,
    Rcpp::NumericVector totmean_yZ );

Rcpp::List mla2_postproc( int N, int NX, int NZ, Rcpp::NumericVector y,
    Rcpp::NumericMatrix X, Rcpp::NumericMatrix Z,
    Rcpp::NumericMatrix idcluster_table,
    Rcpp::NumericVector wgttot, arma::mat theta, arma::mat Tmat,
    arma::mat sig2, arma::mat Var_theta_rj, arma::mat  AfjArj,
    Rcpp::NumericVector wgtlev2, arma::mat Xt2);

Rcpp::List create_dummies_mla2( int GG, Rcpp::NumericVector group,
    Rcpp::NumericMatrix X, Rcpp::NumericMatrix Z, Rcpp::NumericVector y );

Rcpp::NumericMatrix create_idclustertable( Rcpp::NumericVector group,
    Rcpp::NumericVector cluster, int NC);

Rcpp::NumericVector rescale_lev1weights( Rcpp::NumericMatrix idcluster_table,
        Rcpp::NumericVector wgtlev1 );

Rcpp::List mla2_emsteps( Rcpp::NumericMatrix X, Rcpp::NumericMatrix Z,
    Rcpp::NumericMatrix idcluster_table, arma::mat Tmat,
    arma::mat Xa, arma::mat theta, Rcpp::NumericVector y,
    Rcpp::NumericVector wgtlev2, Rcpp::NumericVector wgtlev1,
    arma::mat ArjArj, arma::mat AfjArj, arma::mat sig2,
    Rcpp::NumericVector W1, Rcpp::NumericVector W2,
    arma::mat XtX, arma::mat Xty  );

Rcpp::List mla2_suffstat( Rcpp::NumericMatrix X, Rcpp::NumericMatrix Z,
    Rcpp::NumericVector y,
    Rcpp::NumericVector wgtlev2, Rcpp::NumericVector wgtlev1,
    Rcpp::NumericVector wgttot, Rcpp::NumericMatrix idcluster_table );

Rcpp::List mla2_inits( arma::mat Xa, Rcpp::NumericMatrix X,
    Rcpp::NumericMatrix Z, Rcpp::NumericVector y, int NZ, Rcpp::NumericVector wgtlev1,
    Rcpp::NumericVector wgttot );

Rcpp::List bifie_mla2_estimation( arma::mat theta_init, arma::mat Tmat_init,
        arma::mat sig2_init, int NX, int NZ, int NC, int N,
        Rcpp::NumericMatrix X, Rcpp::NumericMatrix Z,
        Rcpp::NumericVector y, Rcpp::NumericVector wgtlev2,
        Rcpp::NumericVector wgtlev1,
        Rcpp::NumericVector wgttot, Rcpp::NumericMatrix idcluster_table,
        double globconv, int maxiter, Rcpp::NumericMatrix recov_constraint,
        int is_rcov_constraint, int NRC );

Rcpp::List bifie_mla2_estimation_replicates( int N__, int NC__,
    Rcpp::NumericVector wgttot__, Rcpp::NumericMatrix wgtrep__,
    Rcpp::NumericVector wgtlev1__, Rcpp::NumericVector wgtlev2__,
    Rcpp::NumericMatrix idcluster_table2, arma::mat theta0, arma::mat Tmat0,
    arma::mat sig20, int NX, int NZ,  Rcpp::NumericMatrix X__,
    Rcpp::NumericMatrix Z__, Rcpp::NumericVector y__,
    double globconv, int maxiter, int NP,
    Rcpp::NumericMatrix recov_constraint, int is_rcov_constraint, int NRC);

#endif // _BIFIESURVEY_BIFIESURVEY_RCPP_HELPER_H
