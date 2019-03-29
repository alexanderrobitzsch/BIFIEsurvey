//// File Name: bifiesurvey_rcpp_wald_test.cpp
//// File Version: 0.175


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[RcppNOinterfaces(r, cpp)]]

#include <Rcpp.h>
#include <Rmath.h>


using namespace Rcpp;
using namespace arma;

// [include_header_file]
#include "bifiesurvey_rcpp_helper.h"




//**********************************************************
// BIFIE helper function Wald test
Rcpp::List bifiesurvey_rcpp_wald_test_vcov( int VV, Rcpp::NumericVector Ccols,
    Rcpp::NumericMatrix parsM, Rcpp::NumericMatrix parsrepM, int ii, int RR,
    Rcpp::NumericVector fayfac, arma::mat ACdes, arma::colvec Ardes )
{
    double f1=0;
    //*** calculate covariance matrix for imputation ii
    arma::mat var_w = arma::zeros(VV,VV);
    for (int vv1=0;vv1<VV;vv1++){
        for (int vv2=vv1;vv2<VV;vv2++){
            for (int rr=0;rr<RR;rr++){
                f1 = bifiesurvey_rcpp_extract_fayfac( fayfac, rr);
                var_w(vv1,vv2) += f1 * ( parsrepM( Ccols[vv1], rr+ii*RR ) - parsM( Ccols[vv1], ii ) )*
                            ( parsrepM( Ccols[vv2], rr+ii*RR ) - parsM( Ccols[vv2], ii ) );
            }
            var_w(vv2,vv1) = var_w(vv1,vv2);
        }
    }

    //***  compute chi square Wald test statistic
    // compute covariance matrix of hypothesis
    arma::mat var_hyp = arma::mat( ACdes * var_w * arma::trans(ACdes) );

    // compute inverse of variance matrix of hypothesis
    arma::mat var_hypinv = arma::pinv(var_hyp);

    // parameter vector
    arma::colvec parm_vec = arma::zeros(VV,1);
    for (int vv=0;vv<VV;vv++){
        parm_vec(vv,0) = parsM( Ccols[vv], ii );
    }
    // hypothesis statistic
    arma::mat hyp_stat = arma::mat( ACdes * parm_vec - Ardes );
    arma::mat chi2 = arma::mat( arma::trans( hyp_stat ) * var_hypinv * hyp_stat );

    // int df=ACdes.n_rows;
    return Rcpp::List::create(
            Rcpp::Named("chi2") = chi2,
            Rcpp::Named("var_w") = var_w,
            Rcpp::Named("var_hyp") = var_hyp,
            Rcpp::Named("hyp_stat") = hyp_stat
        );
}
//**********************************************************

//*************************************************************************
//  bifiesurvey_rcpp_wald_test
// [[Rcpp::export]]
Rcpp::List bifiesurvey_rcpp_wald_test( Rcpp::NumericMatrix parsM, Rcpp::NumericMatrix parsrepM,
    Rcpp::NumericMatrix Cdes, Rcpp::NumericVector rdes, Rcpp::NumericVector Ccols,
    Rcpp::NumericVector fayfac )
{
    // number of involved variables in the test
    double eps=1e-10;
    int VV = Ccols.size();
    int Nimp = parsM.ncol();
    double Nimp2 = Nimp;
    Nimp2 = Nimp2 + eps;
    int RR = parsrepM.ncol() / Nimp;
    int df = Cdes.nrow();
    Rcpp::NumericMatrix chi2M(Nimp,2);
    arma::mat var_w = arma::zeros(VV,VV);
    arma::mat var_b = arma::zeros(VV,VV);
    Rcpp::NumericVector parsM_sel(VV);
    arma::mat var_wM = arma::zeros(df,df*Nimp);
    Rcpp::NumericMatrix hyp_statM(df,Nimp);

    // arma form of the design matrix Cdes
    arma::mat ACdes = arma::zeros(df,VV);
    for (int dd=0;dd<df;dd++){
        for (int vv=0;vv<VV;vv++){
            ACdes(dd,vv) = Cdes(dd, Ccols[vv] );
        }
    }

    // arma colvec for design vector
    arma::colvec Ardes = arma::zeros(df,1);
    for (int dd=0;dd<df;dd++){
        Ardes(dd,0) = rdes[dd];
    }
    double tmp1=0;
    double tmp2=0;

    for ( int ii=0; ii < Nimp; ii++){  // begin ii
        Rcpp::List res1 = bifiesurvey_rcpp_wald_test_vcov(  VV,  Ccols, parsM,
                    parsrepM, ii, RR, fayfac, ACdes,  Ardes );
        Rcpp::NumericMatrix chi2a = res1["chi2"];
        Rcpp::NumericMatrix var_hyp1 = res1["var_hyp"];
        Rcpp::NumericMatrix hyp_stat1 = res1["hyp_stat"];

        for (int dd=0;dd<df;dd++){
            hyp_statM(dd,ii) = hyp_stat1(dd,0);
            for (int ee=0;ee<df;ee++){
                var_wM(dd,ee+ii*df) = var_hyp1(dd,ee);
            }
        }
        chi2M(ii,0) = chi2a(0,0);
        tmp1 += chi2M(ii,0);
        chi2M(ii,1) = std::sqrt( chi2M(ii,0) );
        tmp2 += chi2M(ii,1);
        Rcpp::NumericMatrix var_w1 = res1["var_w"];
        for (int vv1=0;vv1<VV;vv1++){
            for (int vv2=0;vv2<VV;vv2++){
                var_w(vv1,vv2) += var_w1(vv1,vv2);
            }
        }
    }  // end ii

    ///-------- D2 statistic ----------------
    // calculate ARIV
    double ariv = tmp1 - Nimp * std::pow( tmp2 / Nimp, 2.0 );
    ariv = ariv / ( Nimp - 1 + eps ) * ( 1 + 1 / Nimp2 );

    // calculate D2 statistic
    double D2 = tmp1 / (Nimp*df) - (Nimp2+1)/(Nimp2-1+eps) * ariv;
    D2 = D2 / ( 1 + ariv );

    // calculate degrees of freedom
    double df2 = df;
    double nu3 = 1000;
    if ( Nimp > 1 ){
        nu3 = std::pow( df2, -3.0/Nimp2) * ( Nimp2 - 1 ) *
                std::pow( 1 + 1 / ( ariv + eps ), 2.0 );
    }
    nu3 = bifiesurvey_rcpp_squeeze(nu3, 1, 1000);
    double p_D2 = ::Rf_pf( D2, df, nu3, FALSE, FALSE );

    ///-------- D1 statistic ----------------
    // calculate covariance matrices
    for (int vv1=0;vv1<VV;vv1++){
        for (int vv2=0;vv2<VV;vv2++){
            var_w(vv1,vv2) = var_w(vv1,vv2) / Nimp2;
        }
    }

    // means of all parameters
    for (int vv=0;vv<VV;vv++){
        for (int ii=0;ii<Nimp;ii++){
            parsM_sel[vv] += parsM( Ccols[vv], ii );
        }
        parsM_sel[vv] = parsM_sel[vv] / Nimp2;
    }

    // between matrix
    for (int vv1=0;vv1<VV;vv1++){
        for (int vv2=0;vv2<VV;vv2++){
            for ( int ii=0; ii<Nimp; ii++){
                var_b(vv1,vv2) += ( parsM( Ccols[vv1], ii ) - parsM_sel[vv1] ) *
                                ( parsM( Ccols[vv2], ii ) - parsM_sel[vv2] );
            }
            var_b(vv1,vv2) = 1 / ( Nimp2 - 1 ) * var_b( vv1, vv2 );
        }
    }

    arma::mat ariv_D1a = arma::mat( var_b * arma::pinv(var_w) );
    double ariv_D1 = bifiesurvey_rcpp_arma_trace(ariv_D1a);

    ariv_D1 = ariv_D1 * ( 1 + 1 / Nimp2 ) / df;
    arma::mat var_t1 = arma::mat( (1+ariv_D1) * var_w );

    // hypothesis matrix
    arma::mat var_hyp = arma::mat( ACdes * var_t1 * arma::trans(ACdes) );
    // compute inverse of variance matrix of hypothesis
    arma::mat var_hypinv = arma::pinv(var_hyp);
    // parameter vector
    arma::colvec parm_vec= arma::zeros(VV,1);
    for (int vv=0;vv<VV;vv++){
        parm_vec(vv,0) = parsM_sel[vv];
    }
    // hypothesis statistic
    arma::mat hyp_stat = arma::mat( ACdes * parm_vec - Ardes );
    arma::mat D1 = arma::mat( arma::trans(hyp_stat) * var_hypinv * hyp_stat ) / df;

    // calculate nu2
    double nu2 = 1 + ( 1 - 2 / ( df * Nimp2 - df ) / ariv_D1 );
    nu2 = 4 + ( df * Nimp2 - df - 4 ) * nu2 * nu2;
    if ((Nimp2 < 2) | (nu2 > 1000)){ nu2 = 1000; }
    double D1_stat = D1(0,0);
    double p_D1 = ::Rf_pf( D1_stat, df, nu2, FALSE, FALSE );

    //*************************************************
    // OUTPUT
    return Rcpp::List::create(
            Rcpp::Named("chi2M") = chi2M,
            Rcpp::Named("ariv") = ariv,
            Rcpp::Named("D2") = D2,
            Rcpp::Named("df") = df,
            Rcpp::Named("nu2") = nu2,
            Rcpp::Named("nu3")= nu3,
            Rcpp::Named("p_D1") = p_D1,
            Rcpp::Named("p_D2") = p_D2,
            Rcpp::Named("Nimp") = Nimp,
            Rcpp::Named("RR") = RR,
            Rcpp::Named("fayfac") = fayfac,
            Rcpp::Named("var_w") = var_w,
            Rcpp::Named("var_b") = var_b,
            Rcpp::Named("D1") = D1,
            Rcpp::Named("Ccols") = Ccols,
            Rcpp::Named("parsM_sel") = parsM_sel,
            Rcpp::Named("var_wM") = var_wM,
            Rcpp::Named("hyp_statM") = hyp_statM
        );
}
//*************************************************************************


// Rcpp::Rcout << "ariv_D1 = " <<  ariv_D1 << std::endl;
