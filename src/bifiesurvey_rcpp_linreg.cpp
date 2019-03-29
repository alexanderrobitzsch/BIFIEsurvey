//// File Name: bifiesurvey_rcpp_linreg.cpp
//// File Version: 0.21


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[RcppNOinterfaces(r, cpp)]]

#include <Rcpp.h>
#include <Rmath.h>


using namespace Rcpp;
using namespace arma;

// [include_header_file]
#include "bifiesurvey_rcpp_helper.h"



//*************************************************************************
// linear regression
Rcpp::List bifiesurvey_rcpp_linreg_compute( Rcpp::NumericMatrix dat1,
    Rcpp::NumericVector group_values, Rcpp::NumericVector dep_index,
    Rcpp::NumericVector pre_index, Rcpp::NumericMatrix wgt,
    Rcpp::NumericVector group_index1 )
{

    int WW = wgt.ncol();
    int N = wgt.nrow();
    int VV = pre_index.size();
    int group_index = group_index1[0];
    int GG=group_values.size();
    int VV2=(2*VV+2)*GG;

    Rcpp::NumericMatrix sumwgt1(GG,WW);
    Rcpp::NumericVector ncases(GG);
    Rcpp::NumericMatrix regr_coef(VV2,WW);
    Rcpp::NumericMatrix indcases(N,GG);

    //---- extract usable cases
    for (int nn=0;nn<N;nn++){ // beg nn
        for (int gg=0; gg < GG; gg++ ){  // beg gg
            if ( dat1(nn,group_index) == group_values[gg] ){
                indcases(nn,gg)=1;
            }
            if ( R_IsNA( dat1(nn, dep_index[0] ) ) ){ // beg R_IsNA dep var
                indcases(nn,gg)=0;
            } // end R_IsNA dep var
            for (int vv=0;vv<VV;vv++){ // beg vv  independent vars
                if ( R_IsNA( dat1(nn, pre_index[vv] ) ) ){
                    indcases(nn,gg)=0;
                }
            } // end vv  independent vars
            if ( indcases(nn,gg) == 1 ){
                ncases[gg] ++;
                for (int ww=0;ww<WW; ww++){
                    sumwgt1(gg,ww) += wgt(nn,ww);
                }
                break;
            }
        }   // end gg
    }   // end nn


    double sig2=0;
    double sig2a=0;

    // start computations over groups
    for (int gg=0;gg<GG;gg++){
        int igg=0;
        // create design matrices
        int ngg = ncases[gg];
        arma::mat X0=arma::zeros( ngg, VV );
        arma::mat X=arma::zeros( ngg, VV );
        arma::colvec y0=arma::zeros( ngg );
        arma::colvec y=arma::zeros( ngg );
        arma::mat wgt0=arma::zeros( ngg, WW);
        Rcpp::NumericMatrix M_pre(VV,WW);
        Rcpp::NumericMatrix SD_pre(VV,WW);
        Rcpp::NumericMatrix M_dep(1,WW);
        Rcpp::NumericMatrix SD_dep(1,WW);

        //***** define input matrices
        for (int nn=0; nn<N; nn++){  // beg nn
            if ( indcases(nn,gg)==1){  // beg indcases gg
                y0(igg,0) = dat1(nn, dep_index[0] );
                for (int vv=0;vv<VV;vv++){  // beg vv
                    X0(igg,vv) = dat1(nn, pre_index[vv] );
                }  // end vv
                igg ++;
            } // end if indcases gg
        } // end nn

        //**** define used matrices
        double wtmp=0;
        for ( int ww=0; ww<WW; ww++){ // beg ww
            igg=0;
            for (int nn=0; nn<N; nn++){  // beg nn
                if ( indcases(nn,gg)==1){  // beg indcases gg
                    wgt0(igg,ww) = wgt(nn,ww);
                    wtmp = std::sqrt( wgt(nn,ww));
                    y(igg,0) = y0(igg,0)*wtmp;
                    M_dep(0,ww) += y0(igg,0) * wgt(nn,ww);
                    SD_dep(0,ww) += y0(igg,0) * y0(igg,0) * wgt(nn,ww);
                    for (int vv=0; vv<VV; vv++){  // beg vv
                        X(igg,vv) = X0(igg, vv )*wtmp;
                        M_pre(vv,ww) += X0(igg,vv) * wgt(nn,ww);
                        SD_pre(vv,ww) += X0(igg,vv) * X0(igg,vv) * wgt(nn,ww);
                    }  // end vv
                    igg ++;
                } // end if indcases gg
            } // end nn

            //*** fit linear model
            arma::colvec coef = arma::solve(X, y);      // fit model y ~ X
            // arma::colvec resid = y - X*coef;            // residuals
            // not that the weights are already included in the residual calculation
            arma::colvec resid(ngg);
            for (int hh=0; hh<ngg; hh++){
                resid(hh,0) = y0(hh,0);
                for (int vv=0; vv<VV; vv++){
                    resid(hh,0) = resid(hh,0) - X0(hh,vv) * coef(vv,0);
                }
            }
            // sig2 = arma::as_scalar( arma::trans(resid)*resid );
            sig2=0;
            for (int hh=0; hh<ngg; hh++){
                sig2 += std::pow( resid(hh,0), 2.0) * wgt0(hh,ww);
            }
            double sggww = sumwgt1(gg,ww);
            sig2a = sig2 / ( sggww - VV );

            // collect all regression coefficients
            // unstandardized coefficients
            for (int vv=0; vv<VV; vv++){
                regr_coef( vv + gg*(2*VV+2), ww ) = coef(vv,0);
            }
            regr_coef( VV + gg*(2*VV+2), ww ) = std::sqrt( sig2a );  // sigma
            // compute R^2
            M_dep(0,ww) = M_dep(0,ww) / sggww;
            double sig3 = SD_dep(0,ww) - sggww * std::pow( M_dep(0,ww), 2.0);
            SD_dep(0,ww) = std::sqrt( sig3 / ( sggww - 1 ) );
            regr_coef( VV+1 + gg*(2*VV+2), ww ) = 1 - sig2 / sig3;

            // compute standardized coefficients
            for (int vv=0; vv<VV; vv++){
                M_pre(vv,ww) = M_pre(vv,ww) / sggww;
                SD_pre(vv,ww) = SD_pre(vv,ww) - sggww * std::pow( M_pre(vv,ww), 2.0);
                SD_pre(vv,ww) = std::sqrt( SD_pre(vv,ww) / ( sggww - 1 ) );
                regr_coef( VV+2+vv + gg*(2*VV+2), ww ) = coef(vv,0) / SD_dep(0,ww) * SD_pre(vv,ww);
            }
        } // end ww
    } // end gg

    return Rcpp::List::create(
            Rcpp::Named("ncases") = ncases,
            Rcpp::Named("sumwgt1") = sumwgt1,
            Rcpp::Named("regr_coef") = regr_coef
    );
}
//*************************************************************************


//*************************************************************************
//  bifiesurvey_rcpp_linreg
// [[Rcpp::export]]
Rcpp::List bifiesurvey_rcpp_linreg( Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1,
        Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector dep_index,
        Rcpp::NumericVector pre_index, Rcpp::NumericVector fayfac, Rcpp::NumericVector NI,
        Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values )
{
    int Nimp = NI[0];
    int N = wgt1.nrow();
    int VV = pre_index.size();
    int NV = datalist.ncol();
    int GG=group_values.size();
    Rcpp::NumericMatrix dat1(N,NV);
    int VV2=(2*VV+2)*GG;
    Rcpp::NumericMatrix ncasesM(GG,Nimp);
    Rcpp::NumericMatrix sumwgtM(GG,Nimp);
    Rcpp::NumericMatrix regrcoefM(VV2,Nimp);
    Rcpp::NumericMatrix regrcoef_varM(VV2,Nimp);
    int WW = wgtrep.ncol();
    Rcpp::NumericMatrix regrcoefrepM(VV2,Nimp*WW);
    Rcpp::Rcout << "|";

    //****** loop over imputations
    for (int ii = 0; ii < Nimp; ii++){

        // extract dataset
        dat1 = datalist( Rcpp::Range( ii*N+0, ii*N+ (N-1) ), Rcpp::Range(0,NV-1) );

        //** apply linear regression to one dataset;
        Rcpp::List res1 = bifiesurvey_rcpp_linreg_compute( dat1, group_values,  dep_index,
                                pre_index, wgt1, group_index1 );
        Rcpp::NumericVector ncases = res1["ncases"];
        Rcpp::NumericVector sumwgt0 = matr2vec(res1["sumwgt1"]);
        Rcpp::NumericVector regrcoef0 = matr2vec(res1["regr_coef"]);

        //*** apply linear regression to replicated datasets
        Rcpp::List res2 = bifiesurvey_rcpp_linreg_compute( dat1, group_values,  dep_index, pre_index,
                        wgtrep, group_index1 );
        Rcpp::NumericMatrix sumwgtrep = res2["sumwgt1"];
        Rcpp::NumericMatrix regrcoefrep = res2["regr_coef"];

        // compute standard errors
        Rcpp::NumericVector regrcoef_var = varjack_helper( regrcoef0, regrcoefrep, fayfac );

        for (int gg=0;gg<GG; gg++){
            ncasesM(gg,ii) = ncases[gg];
            sumwgtM(gg,ii) = sumwgt0[gg];
        }
        for (int zz=0;zz<VV2; zz++){
            regrcoefM(zz,ii) = regrcoef0[zz];
            regrcoef_varM(zz,ii) = regrcoef_var[zz];
            for (int ww=0;ww<WW;ww++){
                regrcoefrepM(zz, ww + ii*WW ) = regrcoefrep(zz,ww);
            }
        }
    Rcpp::Rcout << "-" <<  std::flush;
    }  // end ii;  end multiple imputations
    Rcpp::Rcout << "|" << std::endl;

    //*** Rubin inference
    Rcpp::List regrcoefL = rubin_rules_univ( regrcoefM, regrcoef_varM );

    //***** OUTPUT
    return Rcpp::List::create(
            Rcpp::Named("ncasesM") = ncasesM,
            Rcpp::Named("sumwgtM") = sumwgtM,
            Rcpp::Named("regrcoefrepM") = regrcoefrepM,
            Rcpp::Named("regrcoefL") = regrcoefL,
            Rcpp::Named("regrcoefM") = regrcoefM,
            Rcpp::Named("regrcoef_varM") = regrcoef_varM
        );
}
//*************************************************************************
