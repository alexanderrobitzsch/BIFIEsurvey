//// File Name: bifiesurvey_rcpp_logistreg.cpp
//// File Version: 0.23


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
//**** logistic regression
Rcpp::List bifiesurvey_rcpp_logistreg_compute( Rcpp::NumericVector y, Rcpp::NumericMatrix X,
    Rcpp::NumericVector wgt, Rcpp::NumericVector beta0, double eps, int maxiter )
{
    int N=X.nrow();
    int P=X.ncol();
    double t1=0;
    double minval_logit = -15; // minimum value for logit computation

    //*** create matrices in Armadillo
    // design matrix X
    arma::mat Xa0(N,P);
    arma::mat Xa(N,P);
    for (int nn=0;nn<N;nn++){
        for (int pp=0;pp<P;pp++){
            Xa0(nn,pp) = X(nn,pp);
        }
    }
    // outcome matrix y
    arma::colvec ya(N);
    for (int nn=0;nn<N;nn++){
        ya(nn,0) = y[nn];
    }
    // regression coefficients
    arma::colvec beta_old(P);
    arma::colvec beta_new(P);
    for (int pp=0;pp<P;pp++){
        beta_old(pp,0) = beta0[pp];
    }

    // temporary values in iterations
    arma::colvec pred_logit(N);
    arma::colvec ypred(N);
    arma::colvec z(N);
    arma::colvec AM(N);
    arma::colvec wgta(N);
    double pardiff=100;
    int ii=0;

    //--- begin iterations
    while( ( pardiff > eps ) & ( ii < maxiter ) ){

        // calculate predicted logit value and probability
        for( int nn=0; nn<N; nn++){
            pred_logit(nn,0)=0;
            for ( int pp=0; pp <P; pp++){
                pred_logit(nn,0) += Xa0(nn,pp) * beta_old(pp,0);
            }
            if ( pred_logit(nn,0) < minval_logit ){
                pred_logit(nn,0) = minval_logit;
            }
            ypred(nn,0) = 1 / ( 1 + exp( - pred_logit(nn,0) ) );
        }
        // calculate entries for A matrix and outcome z
        for (int nn=0;nn<N;nn++){
            AM(nn,0) = ypred(nn,0) * ( 1 - ypred(nn,0) );
            wgta(nn,0) = std::sqrt( AM(nn,0) * wgt[nn] );
            z(nn,0) = pred_logit(nn,0) + ( ya(nn,0) - ypred(nn,0) )/AM(nn,0);
            z(nn,0) = wgta(nn,0) * z(nn,0);
        }
        for (int nn=0;nn<N;nn++){
            for (int pp=0;pp<P;pp++){
                Xa(nn,pp)=Xa0(nn,pp)*wgta(nn,0);
            }
        }
        // coefficient
        beta_new = arma::solve(Xa, z);      // fit model y ~ X
        // parameter difference
        pardiff=0;
        for (int pp=0;pp<P;pp++){
            t1 = beta_old(pp,0) - beta_new(pp,0);
            if (t1 < 0 ){ t1 = -t1; }
            if (t1 > pardiff){ pardiff = t1; }
        }
        for (int pp=0;pp<P; pp++){
            beta_old(pp,0) = beta_new(pp,0);
        }
        ii ++;
    }  // end iterations

    // R^2 according McKelvey & Zavoina
    double var_logist = 3.289868;
    double var_pred = arma::var(pred_logit);
    double R2 = var_pred / ( var_pred + var_logist );
    Rcpp::NumericVector parm(P+1);
    for (int pp=0;pp<P;pp++){
        parm[pp] = beta_new(pp,0);
    }
    parm[P] = R2;

    //----- OUTPUT
    return Rcpp::List::create(
            Rcpp::Named("pardiff") = pardiff,
            Rcpp::Named("beta") = beta_new,
            Rcpp::Named("parm") = parm,
            Rcpp::Named("iter") = ii
        );
}
//**********************************************************


//*************************************************************************
//  bifiesurvey_rcpp_logistreg
// [[Rcpp::export]]
Rcpp::List bifiesurvey_rcpp_logistreg( Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1,
    Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector dep_index, Rcpp::NumericVector pre_index,
    Rcpp::NumericVector fayfac, Rcpp::NumericVector NI, Rcpp::NumericVector group_index1,
    Rcpp::NumericVector group_values, double eps, int maxiter )
{
    int Nimp = NI[0];
    int RR = wgtrep.ncol();
    int N = wgt1.nrow();
    int VV = pre_index.size();
    int NV = datalist.ncol();
    int group_index = group_index1[0];
    int GG=group_values.size();
    Rcpp::NumericMatrix dat1(N,NV);
    Rcpp::NumericVector tempvec(N);
    int uu=0;
    int VV2=(VV+1)*GG;
    Rcpp::NumericMatrix ncasesM(GG,Nimp);
    Rcpp::NumericMatrix sumwgtM(GG,Nimp);
    Rcpp::NumericMatrix regrcoefM(VV2,Nimp);
    Rcpp::NumericMatrix regrcoef_varM(VV2,Nimp);
    int WW = wgtrep.ncol();
    Rcpp::NumericMatrix regrcoefrepM(VV2,Nimp*WW);
    Rcpp::NumericMatrix tempcoefrepM(VV2,WW);
    Rcpp::NumericVector regrcoef0(VV2);
    Rcpp::NumericVector beta0(VV);
    Rcpp::Rcout << "|";

    // loop over imputed datasets
    for ( int ii=0; ii < Nimp; ii++ ){
        // extract dataset
        dat1 = datalist( Rcpp::Range( ii*N+0, ii*N+(N-1) ), Rcpp::Range(0,NV-1) );

        // loop over group values gg
        int ind=1;
        for (int gg=0; gg < GG; gg ++ ){
            uu=0;
            for (int nn=0;nn<N;nn++){
                ind = 1;
                if ( dat1(nn,group_index) == group_values[gg] ){
                    if ( R_IsNA( dat1(nn,dep_index[0] ) ) ){
                        ind = 0;
                    }  // end NA dep
                    for (int vv=0;vv<VV;vv++){
                        if ( R_IsNA( dat1(nn,pre_index[vv] ) ) ){
                            ind = 0;
                        }  // end NA dep
                    }
                    if ( ind > 0 ){
                        ncasesM(gg,ii) ++;
                        sumwgtM(gg,ii) += wgt1[nn];
                        tempvec[uu] = nn;
                        uu ++;
                    }
                } // end group_val == gg
            }

            // create datasets for logistic regression
            int ngg=ncasesM(gg,ii);
            Rcpp::NumericMatrix Xt(ngg,VV);
            Rcpp::NumericVector yt(ngg);
            Rcpp::NumericVector wgtt(ngg);
            Rcpp::NumericMatrix wgtrept(ngg,RR);
            Rcpp::NumericVector wgttemp(ngg);

            for (int tt=0;tt<ngg;tt++){
                yt[tt] = dat1( tempvec[tt], dep_index[0] );
                wgtt[tt] = wgt1[ tempvec[tt] ];
                for (int rr=0;rr<RR;rr++){
                    wgtrept(tt,rr) = wgtrep( tempvec[tt],rr);
                }
                for (int vv=0;vv<VV;vv++){
                    Xt(tt,vv)=dat1( tempvec[tt], pre_index[vv] );
                }
            } // end cases tt

            // logistic regression original dataset
            Rcpp::List res1 = bifiesurvey_rcpp_logistreg_compute( yt, Xt, wgtt,
                                    beta0, eps, maxiter );
            Rcpp::NumericVector tempcoef = res1["beta"];
            Rcpp::NumericVector tempparm = res1["parm"];
            for (int vv=0;vv<VV+1;vv++){
                regrcoefM(vv+gg*VV,ii) = tempparm[vv];
            }

            Rcpp::List res2;
            // replicated datasets
            for (int rr=0;rr<RR;rr++){
                for (int tt=0;tt<ngg;tt++){
                    wgttemp[tt] = wgtrept(tt,rr);
                }
                res2 = bifiesurvey_rcpp_logistreg_compute( yt, Xt, wgttemp,
                            tempcoef, eps, maxiter );
                Rcpp::NumericVector tempcoef2=res2["parm"];
                for (int vv=0;vv<VV+1;vv++){
                    tempcoefrepM(vv+gg*VV,rr) = tempcoef2[vv];
                }
            } // end rr
        } //----- end gg

        for (int zz=0;zz<VV2;zz++){
            regrcoef0[zz] = regrcoefM(zz,ii);
        }

        // compute standard errors
        Rcpp::NumericVector regrcoef_var = varjack_helper( regrcoef0,
                        tempcoefrepM, fayfac );

        for (int zz=0;zz<VV2; zz++){
            regrcoef_varM(zz,ii) = regrcoef_var[zz];
            for (int ww=0;ww<WW;ww++){
                regrcoefrepM(zz, ww + ii*WW ) = tempcoefrepM(zz,ww);
            }
        }
        Rcpp::Rcout << "-" <<  std::flush;
    }  // end ii;  end multiple imputations

    Rcpp::Rcout << "|" << std::endl;
    ///*** Rubin inference
    Rcpp::List regrcoefL = rubin_rules_univ( regrcoefM, regrcoef_varM );

    //----- OUTPUT
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
