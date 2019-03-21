//// File Name: bifiesurvey_rcpp_pathmodel.cpp
//// File Version: 0.14


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[RcppNOinterfaces(r, cpp)]]

#include <Rcpp.h>
#include <Rmath.h>


using namespace Rcpp;
using namespace arma;

// [include_header_file]
#include "bifiesurvey_rcpp_helper.h"


//***************************************
// path model
Rcpp::List bifie_rcpp_pathmodel_compute( Rcpp::NumericMatrix dat1,
    Rcpp::NumericMatrix wgt, Rcpp::NumericVector group_values,
    Rcpp::NumericVector group_index1, Rcpp::NumericVector vars_index,
    Rcpp::NumericMatrix R, int NP0, Rcpp::NumericMatrix coeff_index,
    Rcpp::NumericMatrix E, Rcpp::NumericMatrix L, int NL,
    Rcpp::NumericVector L_row_index, Rcpp::NumericVector R_row_index,
    Rcpp::NumericVector unreliability )
{

    int N = dat1.nrow();
    int GG = group_values.size();
    int RR = wgt.ncol();
    int VV = vars_index.size();

    int group_index = group_index1[0];
    Rcpp::NumericMatrix indcases(N,GG);
    Rcpp::NumericVector ncases(GG);
    Rcpp::NumericMatrix sumwgt1(GG,RR);
    int WW = wgt.ncol();
    int NR = R.nrow();

    // parameter vector
    int NP = 2*NP0 + 2*NR;
    Rcpp::NumericMatrix parsM( GG*NP, WW );
    int Npow = coeff_index.ncol();
    int Nrci = coeff_index.nrow();

    //***** extract usable cases
    for (int nn=0;nn<N;nn++){ // beg nn
        for (int gg=0; gg < GG; gg++ ){  // beg gg
            if ( dat1(nn,group_index) == group_values[gg] ){
                indcases(nn,gg)=1;
            }
            for (int vv=0;vv<VV;vv++){ // beg vv  independent vars
                if ( R_IsNA( dat1(nn, vars_index[vv] ) ) ){
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


    // create matrices for group gg
    int igg=0;
    // means, SDs and covariances
    Rcpp::NumericMatrix means_gg(VV,GG);
    Rcpp::NumericMatrix sds_gg(VV,GG);
    Rcpp::NumericVector weights_gg(GG);
    Rcpp::NumericMatrix covs_gg(VV,VV*GG);
    Rcpp::NumericMatrix XtX_gg(VV,VV*GG);

    Rcpp::NumericVector rsquared(NR);
    Rcpp::NumericVector residvar(NR);
    Rcpp::NumericVector preds_rr(VV);

    for ( int gg=0; gg<GG;gg++){
        igg=0;
        int ngg = ncases[gg];
        arma::mat X=arma::zeros(ngg,VV);
        Rcpp::NumericMatrix wgg(ngg,WW);
        Rcpp::NumericMatrix E_gg(VV,VV);
        for (int vv=0;vv<VV;vv++){
            E_gg(_,vv) = E(_,vv);
        }

        for (int nn=0;nn<N;nn++){
            if ( indcases(nn,gg) == 1){
                for (int vv=0;vv<VV;vv++){
                    if ( ! R_IsNA(vars_index[vv] ) ){
                        X(igg,vv) = dat1( nn, vars_index[vv] );
                    }
                }
                for (int ww=0;ww<WW;ww++){
                    wgg(igg,ww) = wgt(nn,ww);
                }
                igg ++;
            }
        }

        // create derived variables
        if (NL>0){
            for (int jj=0;jj<NL;jj++){
                for (int nn=0;nn<ngg;nn++){
                    for (int vv=0;vv<VV;vv++){
                        if ( L(jj,vv) != 0 ){
                            X( nn, L_row_index[jj] ) += X( nn, vv ) * L( jj, vv );
                        }
                    }
                }
            }
        }

        // calculate means, variances and covaiances
        double cbar=0;
        double vbar=0;
        double alpha=0;
        double NItems=0;
        int lli=0;
        for (int ww=0;ww<WW;ww++){
            for (int vv=0;vv<VV;vv++){
                sds_gg(vv,gg)=0;
                means_gg(vv,gg)=0;
                for (int xx=0;xx<VV;xx++){
                    covs_gg(vv,xx+gg*VV)=0;
                }
            }

            weights_gg[gg] = 0;
            arma::mat B_est_gg =arma::zeros(VV,VV);
            arma::mat B_est_stand_gg =arma::zeros(VV,VV);

            for (int nn=0;nn<ngg;nn++){
                weights_gg[gg] += wgg(nn,ww);
                for (int vv=0;vv<VV;vv++){
                    means_gg(vv,gg) += wgg(nn,ww) * X(nn,vv);
                    sds_gg(vv,gg) += wgg(nn,ww) * X(nn,vv) * X(nn,vv);
                    for (int xx=vv+1;xx<VV;xx++){
                        covs_gg( vv, xx+gg*VV) += wgg(nn,ww) * X(nn,vv) * X(nn,xx);
                    }
                }
            }
            for (int vv=0;vv<VV;vv++){
                means_gg(vv,gg) = means_gg( vv,gg) / weights_gg[gg];
                sds_gg(vv,gg) = std::sqrt( ( sds_gg(vv,gg) - weights_gg[gg] * std::pow( means_gg( vv, gg), 2.0 ) )/
                        weights_gg[gg] );
                covs_gg(vv,vv+gg*VV) = std::pow( sds_gg(vv,gg), 2.0 );
            }

            for (int vv=0;vv<VV;vv++){
                for (int xx=vv+1;xx<VV;xx++){
                    covs_gg( vv, xx+gg*VV) = covs_gg(vv,xx+gg*VV) - weights_gg[gg] * means_gg(vv,gg) * means_gg(xx,gg);
                    covs_gg( vv, xx+gg*VV) = covs_gg(vv,xx+gg*VV) / weights_gg[gg];
                    covs_gg( xx, vv+gg*VV) = covs_gg(vv,xx+gg*VV);
                }
            }

            // calculate Cronbach's alpha for latent variables and adjust
            // measurement error variance
            for (int jj=0;jj<NL;jj++){
                NItems=0;
                for (int vv=0;vv<VV;vv++){
                    vbar += std::sqrt( L(jj,vv) * L(jj,vv) ) * covs_gg( vv, vv+gg*VV);
                    if ( L(jj,vv) != 0 ){
                        NItems = NItems + 1;
                    }
                    for (int xx=0;xx<VV;xx++){
                        if ( xx != vv ){
                            cbar += L(jj,vv) * L(jj,xx) * covs_gg( vv, xx+gg*VV);
                        }
                    }
                }
                vbar = vbar / NItems;
                cbar = cbar / ( NItems*(NItems-1) );
                alpha = NItems * cbar / ( vbar + ( NItems - 1 ) * cbar );
                lli = L_row_index[jj];
                E_gg( lli, lli ) = ( 1 - alpha ) * covs_gg( lli, lli +gg*VV);
            }

            // correction for unreliability of variances
            for (int vv=0;vv<VV;vv++){
                if ( unreliability[vv] > 0 ){
                    E_gg( vv, vv ) = covs_gg( vv, vv + gg*VV ) * unreliability[vv];
                }
            }

            // adjustment of covariance matrix for measurement error covariance
            for (int vv=0;vv<VV;vv++){
                for (int xx=0;xx<VV;xx++){
                    XtX_gg(vv,xx+gg*VV) = covs_gg( vv, xx+gg*VV) - E_gg(vv,xx);
                }
            }

            // calculate regression coefficients
            int nv_rr = 0;

            for (int rr=0;rr<NR;rr++){
                // select predictors
                nv_rr=0;
                for (int vv=0;vv<VV;vv++){
                    if ( ( R(rr,vv) != 0 ) ){
                        preds_rr[nv_rr] = vv;
                        nv_rr = nv_rr+1;
                    }
                }
                arma::mat sigma_XX=arma::zeros(nv_rr,nv_rr);
                arma::mat sigma_XY=arma::zeros(nv_rr,1);
                arma::mat regr_coef_stand=arma::zeros(nv_rr,1);
                for (int vv=0;vv < nv_rr; vv++){
                    for (int xx=0;xx < nv_rr; xx++){
                        sigma_XX( vv, xx ) = XtX_gg( preds_rr[vv], preds_rr[xx]+gg*VV);
                    }
                    sigma_XY(vv,0) = XtX_gg( preds_rr[vv], R_row_index[rr]+gg*VV);
                }

                arma::mat regr_coef = arma::mat( arma::pinv( sigma_XX ) * sigma_XY );
                for (int vv=0; vv < nv_rr; vv ++ ){
                    regr_coef_stand(vv,0) = regr_coef(vv,0) * std::sqrt( sigma_XX(vv,vv) ) /
                            std::sqrt( XtX_gg( R_row_index[rr], R_row_index[rr]+gg*VV) );
                    B_est_gg( R_row_index[rr], preds_rr[vv] ) = regr_coef(vv,0);
                    B_est_stand_gg( R_row_index[rr], preds_rr[vv] ) = regr_coef_stand(vv,0);
                }
                // calculate R-squared
                arma::mat var_expl = arma::mat( arma::trans( regr_coef ) * sigma_XX * regr_coef );
                rsquared[rr] = var_expl(0,0) / XtX_gg( R_row_index[rr], R_row_index[rr]+gg*VV);
                residvar[rr] = XtX_gg( R_row_index[rr], R_row_index[rr]+gg*VV) * ( 1- rsquared[rr] );
            }

            //--- calculate parameter list
            // extract coefficients
            for (int zz=0;zz<Nrci; zz++){ // beg zz
                parsM( zz + gg*NP, ww ) = 1;
                parsM( zz + NP0 + gg*NP, ww ) = 1;
                for (int hh=2;hh<Npow-1;hh++){
                    if ( ! ( R_IsNA( coeff_index(zz,hh+1) ) ) ){
                        parsM( zz + gg*NP, ww ) = parsM( zz + gg*NP, ww ) *
                        B_est_gg( coeff_index(zz,hh+1)-1, coeff_index(zz,hh)-1 );
                        parsM( zz + NP0+ gg*NP, ww ) = parsM( zz + NP0 + gg*NP, ww ) *
                        B_est_stand_gg( coeff_index(zz,hh+1)-1, coeff_index(zz,hh)-1 );
                    }
                }
                if ( ! ( R_IsNA( coeff_index(zz, 0 ) ) ) ){
                    parsM( coeff_index(zz, 0 ) - 1 + gg*NP, ww ) += parsM( zz + gg*NP, ww );
                    parsM( coeff_index(zz, 0 ) - 1 + NP0 + gg*NP, ww ) += parsM( zz + NP0+ gg*NP, ww );
                }
                if ( ! ( R_IsNA( coeff_index(zz, 1 ) ) ) ){
                    parsM( coeff_index(zz, 1 ) - 1 + gg*NP, ww ) += parsM( zz + gg*NP, ww );
                    parsM( coeff_index(zz, 1 ) - 1 + NP0 + gg*NP, ww ) += parsM( zz + NP0+ gg*NP, ww );
                }
            }  // end zz
            // r squared
            for (int rr=0;rr<NR;rr++){
                parsM( rr + 2*NP0 + gg*NP, ww ) = rsquared[rr];
                parsM( rr + NR + 2*NP0 + gg*NP, ww ) = residvar[rr];
            }
        } // end ww
    }  // end gg

    //--- output
    return Rcpp::List::create(
            Rcpp::Named("parsM") = parsM,
            Rcpp::Named("ncases") = ncases,
            Rcpp::Named("sumwgt1") = sumwgt1
        );
}
//**********************************************************


//*************************************************************************
//  bifiesurvey_rcpp_pathmodel
// [[Rcpp::export]]
Rcpp::List bifiesurvey_rcpp_pathmodel( Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1,
        Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector vars_index,
        Rcpp::NumericVector fayfac, Rcpp::NumericVector NI,
        Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values,
        Rcpp::NumericMatrix L, Rcpp::NumericVector L_row_index,
        int NL, Rcpp::NumericMatrix E, Rcpp::NumericMatrix R, Rcpp::NumericVector R_row_index,
        Rcpp::NumericMatrix coeff_index, int NP0, Rcpp::NumericVector unreliability )
{

    int Nimp = NI[0];
    int RR = wgtrep.ncol();
    int N = wgt1.nrow();
    int NV = datalist.ncol();
    int GG=group_values.size();

    // parameter vector
    int NR = R.nrow();
    int NP = 2*NP0 + 2*NR;
    Rcpp::NumericMatrix dat1(N,NV);
    Rcpp::NumericMatrix parsM( NP*GG, Nimp );
    Rcpp::NumericMatrix ncases( GG, 1 );
    Rcpp::NumericMatrix sumwgt( GG, 1 );
    Rcpp::NumericMatrix parsrepM( NP*GG, Nimp*RR );
    Rcpp::NumericMatrix parsVar( NP*GG, Nimp);
    Rcpp::Rcout << "|";

    // loop imputed datasets
    for (int ii=0;ii<Nimp;ii++){
        //--- extract dataset
        dat1 = datalist( Rcpp::Range( ii*N+0, ii*N+ (N-1) ), Rcpp::Range(0,NV-1) );
        //--- computation complete data
        Rcpp::List res21 = bifie_rcpp_pathmodel_compute( dat1, wgt1, group_values, group_index1,
                vars_index, R, NP0, coeff_index, E, L, NL, L_row_index,
                R_row_index, unreliability );
        Rcpp::NumericMatrix pM = res21["parsM"];
        parsM(_,ii) = pM(_,0);
        if (ii==0){
            Rcpp::NumericVector v1 = res21["ncases"];
            ncases(_,0) = v1;
            Rcpp::NumericMatrix v2 = res21["sumwgt1"];
            sumwgt(_,0) = v2(_,0);
        }

        //--- computation replicated data
        Rcpp::List res22 = bifie_rcpp_pathmodel_compute( dat1, wgtrep, group_values, group_index1,
                vars_index, R, NP0, coeff_index, E, L, NL, L_row_index,
                R_row_index, unreliability );
        Rcpp::NumericMatrix pM1 = res22["parsM"];
        for (int rr=0;rr<RR; rr++){
            parsrepM( _, rr + ii*RR ) = pM1( _, rr );
        }

        //--- replication variance
        Rcpp::List res41 = varjack_bias_helper( pM, pM1, fayfac );
        Rcpp::NumericVector pars_var = res41["pars_var"];
        int NP1 = pars_var.size();
        for (int pp=0;pp<NP1;pp++){
            parsVar( pp, ii ) = pars_var[ pp ];
        }
        Rcpp::Rcout << "-" <<  std::flush;
    } // end ii (imputations)

    //*** Rubin inference
    Rcpp::List parsL = rubin_rules_univ( parsM, parsVar );
    Rcpp::Rcout << "|" << std::endl;

    //***** OUTPUT
    return Rcpp::List::create(
            Rcpp::Named("parsL") = parsL,
            Rcpp::Named("parsM") = parsM,
            Rcpp::Named("ncases") = ncases,
            Rcpp::Named("sumwgt") = sumwgt,
            Rcpp::Named("parsrepM") = parsrepM,
            Rcpp::Named("parsVar") = parsVar
        );
}
//*************************************************************************

