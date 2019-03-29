//// File Name: bifiesurvey_rcpp_main.cpp
//// File Version: 7.946


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[RcppNOinterfaces(r, cpp)]]

#include <Rcpp.h>
#include <Rmath.h>


using namespace Rcpp;
using namespace arma;

// [include_header_file]
#include "bifiesurvey_rcpp_helper.h"

// [include_header_file]
#include "bifiesurvey_rcpp_wald_test.h"



//*************************************************************************
//  univar_multiple_V2group
// [[Rcpp::export]]
Rcpp::List univar_multiple_V2group( Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1,
      Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector vars_index, Rcpp::NumericVector fayfac,
      Rcpp::NumericVector NI, Rcpp::NumericVector group_index1,
      Rcpp::NumericVector group_values )
{

     int Nimp = NI[0];
     int RR = wgtrep.ncol();
     int N = wgt1.nrow();
     int VV = vars_index.size();
     int NV = datalist.ncol();
     int GG=group_values.size();
     int VV1=VV*GG;

     Rcpp::NumericMatrix mean1M(VV1,Nimp);
     Rcpp::NumericMatrix sd1M(VV1,Nimp);
     Rcpp::NumericMatrix mean1_varM(VV1,Nimp);
     Rcpp::NumericMatrix sd1_varM(VV1,Nimp);
     Rcpp::NumericMatrix mean1repM(VV1,RR*Nimp);
     Rcpp::NumericMatrix sd1repM(VV1,RR*Nimp);
     Rcpp::NumericMatrix dat1(N,NV);
     Rcpp::NumericVector sumweights(1);
     Rcpp::NumericVector mean1(VV1);
     Rcpp::NumericVector sd1(VV1);
     Rcpp::NumericVector sumwgt1(VV);
     Rcpp::NumericVector ncases1(VV*GG);
     Rcpp::NumericVector mean1_var(VV1);
     Rcpp::NumericVector sd1_var(VV1);
     Rcpp::NumericMatrix sumweightM(GG,Nimp);
     Rcpp::NumericMatrix ncasesM(GG*VV,Nimp);
     Rcpp::NumericMatrix mean1rep(VV1,RR);
     Rcpp::NumericMatrix sd1rep(VV1,RR);
     Rcpp::NumericMatrix sumweightrepM(GG,RR*Nimp);

     Rcpp::Rcout << "|";

     //***********************
     // loop multiply imputed datasets
     for ( int ii=0; ii < Nimp; ii++){

     dat1 = datalist( Rcpp::Range( ii*N+0, ii*N+(N-1) ), Rcpp::Range(0,NV-1) );
     // statistics single data
     Rcpp::List res1 = univar_helper_multiple_V2group( dat1,  wgt1,  vars_index,
         group_index1, group_values );
     mean1 = matr2vec(res1["mean1"]);
     sd1 = matr2vec(res1["sd1"]);
     sumwgt1 = matr2vec(res1["sumwgt1"]);
     ncases1 = matr2vec(res1["ncases1"]);

     // compute statistics for replicate weights;
     Rcpp::List res2 = univar_helper_multiple_V2group( dat1,  wgtrep,  vars_index,
         group_index1, group_values );
     Rcpp::NumericMatrix mean1rep = res2["mean1"];
     Rcpp::NumericMatrix sd1rep = res2["sd1"];
//     Rcpp::NumericMatrix sd1rep = res2["sd1"];
     Rcpp::NumericMatrix sumwgt1a = res2["sumwgt1"];

         // compute standard errors
         mean1_var = varjack_helper( mean1, mean1rep, fayfac );
         sd1_var = varjack_helper( sd1, sd1rep, fayfac );

         // collect all results for one imputed dataset
         for (int vv=0;vv<VV1;vv++){
            mean1M( vv, ii ) = mean1[vv];
            sd1M( vv, ii ) = sd1[vv];
            mean1_varM( vv, ii ) = mean1_var[vv];
            sd1_varM( vv, ii ) = sd1_var[vv];
            for (int rr=0;rr<RR;rr++){
             mean1repM( vv, rr+ii*RR ) = mean1rep(vv,rr);
             sd1repM( vv, rr+ii*RR ) = sd1rep(vv,rr);
                           }
                 }

            for (int gg=0;gg<GG;gg++){
                      sumweightM(gg,ii) = sumwgt1[gg];
                      for (int rr=0;rr<RR;rr++){
                          sumweightrepM( gg, rr+ii*RR) = sumwgt1a(gg,rr);
                      }
                      for (int vv=0;vv<VV;vv++){
                         ncasesM(gg+vv*GG,ii) = ncases1[gg+vv*GG];
                                 }
                                   }
           Rcpp::Rcout << "-" <<  std::flush;
            // << std::endl;
         } // end loop ii | imputed datasets

        Rcpp::Rcout << "|" << std::endl;

     //----
     // inference multiply imputed datasets
     Rcpp::List res3 = rubin_rules_univ( mean1M, mean1_varM );
     mean1=res3["pars"];
     Rcpp::NumericVector mean1_se=res3["pars_se"];
     Rcpp::NumericVector mean1_varWithin=res3["pars_varWithin"];
     Rcpp::NumericVector mean1_varBetween=res3["pars_varBetween"];
     Rcpp::NumericVector mean1_fmi=res3["pars_fmi"];
     res3 = rubin_rules_univ( sd1M, sd1_varM );
     sd1=res3["pars"];
     Rcpp::NumericVector sd1_se=res3["pars_se"];
     Rcpp::NumericVector sd1_varWithin=res3["pars_varWithin"];
     Rcpp::NumericVector sd1_varBetween=res3["pars_varBetween"];
     Rcpp::NumericVector sd1_fmi=res3["pars_fmi"];

     res3 = rubin_rules_univ( ncasesM, ncasesM );
     Rcpp::NumericVector ncases =res3["pars"];

     //*************************************************
     // OUTPUT
     return Rcpp::List::create(
         Rcpp::_["mean1"] = mean1,
         Rcpp::_["mean1_se"] = mean1_se,
         Rcpp::_["mean1_varWithin"] = mean1_varWithin,
         Rcpp::_["mean1_varBetween"] = mean1_varBetween,
         Rcpp::_["mean1_fmi"] = mean1_fmi,
         Rcpp::_["mean1M"] = mean1M,
         Rcpp::_["mean1_varM"] = mean1_varM,
         Rcpp::_["mean1repM"] = mean1repM,
         Rcpp::_["sd1"] = sd1,
         Rcpp::_["sd1_se"] = sd1_se,
         Rcpp::_["sd1_varWithin"] = sd1_varWithin,
         Rcpp::_["sd1_varBetween"] = sd1_varBetween,
         Rcpp::_["sd1_fmi"] = sd1_fmi,
         Rcpp::_["sd1M"] = sd1M,
         Rcpp::_["sd1_varM"] = sd1_varM,
         Rcpp::_["sd1repM"] = sd1repM,
         Rcpp::_["sumweightrepM"] = sumweightrepM,
         Rcpp::_["sumweightM"] = sumweightM,
         Rcpp::_["ncases"] = ncases,
         Rcpp::_["ncasesM"] = ncasesM
         );
}

//*************************************************************************


//*************************************************************************
//  bifie_freq
// [[Rcpp::export]]
Rcpp::List bifie_freq( Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1,
    Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector vars_index,
    Rcpp::NumericVector fayfac, Rcpp::NumericVector NI, Rcpp::NumericVector group_index1,
    Rcpp::NumericVector group_values, Rcpp::NumericMatrix vars_values,
    Rcpp::NumericVector vars_values_numb ){

     int Nimp = NI[0];
     int RR = wgtrep.ncol();
     int N = wgt1.nrow();
     int VV = vars_index.size();
     int NV = datalist.ncol();
     // int group_index = group_index1[0];
     int GG=group_values.size();

     // number of values per variable
     int VV1=0;
     Rcpp::NumericVector vars_values_numb_cumsum(VV+1);
     for (int vv=0;vv<VV; vv++){
        VV1 += vars_values_numb[vv];
        vars_values_numb_cumsum[vv+1] = vars_values_numb_cumsum[vv] + vars_values_numb[vv];
                }
     int VV2=VV1*GG;

     // matrices for output
     int WW = wgt1.ncol();

     Rcpp::NumericMatrix dat1(N,NV);
     Rcpp::NumericMatrix perc1(VV2,WW);
     Rcpp::NumericMatrix perc2(VV2,WW);
     Rcpp::NumericVector sumwgt(VV*GG);
     Rcpp::NumericVector ncases(VV*GG);
     Rcpp::NumericVector ncases1(VV2);
     Rcpp::NumericVector perc1a(VV2);
     Rcpp::NumericVector perc2a(VV2);
     Rcpp::NumericMatrix perc1M(VV2,Nimp);
     Rcpp::NumericMatrix perc2M(VV2,Nimp);
     Rcpp::NumericMatrix perc1_varM(VV2,Nimp);
     Rcpp::NumericMatrix perc2_varM(VV2,Nimp);
     Rcpp::NumericMatrix ncases1M(VV2,Nimp);
     Rcpp::NumericMatrix perc1repM(VV2,RR*Nimp);
     Rcpp::NumericMatrix perc2repM(VV2,RR*Nimp);

     Rcpp::Rcout << "|";

     // loop over imputed datasets
     for (int ii = 0;ii<Nimp; ii++ ){

     dat1 = datalist( Rcpp::Range( ii*N+0, ii*N+ (N-1) ), Rcpp::Range(0,NV-1) );

     //****** statistics single data
     Rcpp::List res1 = bifiehelper_freq( dat1, wgt1, group_index1, group_values,
         vars_values_numb, vars_values, vars_index,  vars_values_numb_cumsum );
     perc1a = matr2vec(res1["perc1"]);
     perc2a = matr2vec(res1["perc2"]);
     ncases = res1["ncases"];
     ncases1 = res1["ncases1"];
     sumwgt = matr2vec(res1["sumwgt"]);

     //****** statistics replicated dataset
     Rcpp::List res2 = bifiehelper_freq( dat1, wgtrep, group_index1, group_values,
         vars_values_numb, vars_values, vars_index,  vars_values_numb_cumsum );

     // compute statistics for replicate weights;
     Rcpp::NumericMatrix perc1rep = res2["perc1"];
     Rcpp::NumericMatrix perc2rep = res2["perc2"];

     // compute standard errors
     Rcpp::NumericVector perc1_var = varjack_helper( perc1a, perc1rep, fayfac );
     Rcpp::NumericVector perc2_var = varjack_helper( perc2a, perc2rep, fayfac );

     for (int zz=0;zz<VV2;zz++){
          perc1M(zz,ii) = perc1a[zz];
          perc2M(zz,ii) = perc2a[zz];
          perc1_varM(zz,ii) = perc1_var[zz];
          perc2_varM(zz,ii) = perc2_var[zz];
          ncases1M(zz,ii) = ncases1[zz];
          for (int rr=0;rr<RR;rr++){
             perc1repM( zz, rr+ii*RR ) = perc1rep(zz,rr);
             perc2repM( zz, rr+ii*RR ) = perc2rep(zz,rr);
           }
     }

     Rcpp::Rcout << "-" <<  std::flush;

     } // end ii

     Rcpp::Rcout << "|" << std::endl;

     ///*** Rubin inference
     Rcpp::List perc1L = rubin_rules_univ( perc1M, perc1_varM );
     Rcpp::List perc2L = rubin_rules_univ( perc2M, perc2_varM );

     // another output list
     Rcpp::List out1 = Rcpp::List::create(
         Rcpp::_["GG"] = GG,
         Rcpp::_["VV2"] = VV2
             );
     //*************************************************
     // OUTPUT
     return Rcpp::List::create(
         Rcpp::_["ncases1M"] = ncases1M,
         Rcpp::_["ncases"] = ncases,
         Rcpp::_["perc1"] = perc1L,
         Rcpp::_["perc1M"] = perc1M,
         Rcpp::_["perc1_varM"] = perc1_varM,
         Rcpp::_["perc1repM"] = perc1repM,
         Rcpp::_["perc2"] = perc2L,
         Rcpp::_["perc2M"] = perc2M,
         Rcpp::_["perc2_varM"] = perc2_varM,
         Rcpp::_["perc2repM"] = perc2repM,
         Rcpp::_["outlist"] = out1
         );
}
//*************************************************************************



//*************************************************************************
//  bifie_correl
// [[Rcpp::export]]
Rcpp::List bifie_correl( Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1,
    Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector vars_index,
    Rcpp::NumericVector fayfac, Rcpp::NumericVector NI, Rcpp::NumericVector group_index1,
    Rcpp::NumericVector group_values ){

     int Nimp = NI[0];
     // int RR = wgtrep.ncol();

     int N = wgt1.nrow();
     int VV = vars_index.size();
     int NV = datalist.ncol();
     // int group_index = group_index1[0];
     int GG=group_values.size();

     Rcpp::NumericMatrix dat1(N,NV);
     int WW = wgtrep.ncol();
     int RR = WW;
     Rcpp::NumericMatrix mean1(VV*GG,WW);
     Rcpp::NumericMatrix sd1(VV*GG,WW);
     Rcpp::NumericVector sumwgt1(GG);
     Rcpp::NumericVector ncases1(GG);
     Rcpp::NumericMatrix ncases1M(GG,Nimp);
     Rcpp::NumericMatrix sumwgt1M(GG,Nimp);
     // create index matrix for covariances and correlations
     int ZZ = VV*(VV-1) / 2 + VV;
     Rcpp::NumericMatrix itempair_index( ZZ, 2 );
     Rcpp::NumericVector cov1(ZZ*GG);
     Rcpp::NumericVector cor1(ZZ*GG);

     int zz=0;
     for (int vv1=0;vv1<VV;vv1++){
       for (int vv2=vv1;vv2<VV;vv2++){
          itempair_index(zz,0) = vv1;
          itempair_index(zz,1) = vv2;
          zz++;
       }
    }

     int VV2 = ZZ*GG;
     Rcpp::NumericMatrix cov1M(VV2,Nimp);
     Rcpp::NumericMatrix cov1_varM(VV2,Nimp);
     Rcpp::NumericMatrix cov1repM(VV2,RR*Nimp);
     Rcpp::NumericMatrix cor1M(VV2,Nimp);
     Rcpp::NumericMatrix cor1_varM(VV2,Nimp);
     Rcpp::NumericMatrix cor1repM(VV2,RR*Nimp);
     Rcpp::Rcout << "|";

     ///***************** loop imputed datasets

     for (int ii = 0; ii < Nimp; ii++){
     // extract dataset
     dat1 = datalist( Rcpp::Range( ii*N+0, ii*N+ (N-1) ), Rcpp::Range(0,NV-1) );

     // perform listwise deletion of cases
     // create indicators of persons
     Rcpp::NumericVector ind_cases(N);
     for ( int nn=0; nn<N; nn++){ // beg nn
       ind_cases[nn] = 1;
       for ( int vv=0; vv<VV; vv++){ // beg vv
           if ( R_IsNA( dat1(nn, vars_index[vv] ) ) ){
                       ind_cases[nn] = 0;
                           }
                    }  // end vv
                  } // end nn

     // statistics one dataset
     Rcpp::List res2 = bifiehelpers_correl( dat1, ind_cases, group_values, group_index1,
          wgt1,  vars_index, itempair_index );
     // extract statistics
     cov1 = matr2vec( res2["cov1"] );
     cor1 = matr2vec( res2["cor1"] );
     sumwgt1 = res2["sumwgt1"];
     ncases1 =  res2["ncases1"];

     // statistics replicated datasets
     Rcpp::List res3 = bifiehelpers_correl( dat1, ind_cases, group_values, group_index1,
          wgtrep,  vars_index, itempair_index );

     // compute statistics for replicate weights;
     Rcpp::NumericMatrix cov1rep = res3["cov1"];
     Rcpp::NumericMatrix cor1rep = res3["cor1"];

     // compute standard errors
     Rcpp::NumericVector cov1_var = varjack_helper( cov1, cov1rep, fayfac );
     Rcpp::NumericVector cor1_var = varjack_helper( cor1, cor1rep, fayfac );

     for (int zz=0;zz<VV2;zz++){
          cov1M(zz,ii) = cov1[zz];
          cov1_varM(zz,ii) = cov1_var[zz];
          cor1M(zz,ii) = cor1[zz];
          cor1_varM(zz,ii) = cor1_var[zz];
          for (int rr=0;rr<RR;rr++){
             cov1repM( zz, rr+ii*RR ) = cov1rep(zz,rr);
             cor1repM( zz, rr+ii*RR ) = cor1rep(zz,rr);
                           }
               }
     for (int gg=0;gg<GG;gg++){
          ncases1M(gg,ii) = ncases1[gg];
          sumwgt1M(gg,ii) = sumwgt1[gg];
                  }

     Rcpp::Rcout << "-" <<  std::flush;

          }  // end ii;  end multiple imputations

     Rcpp::Rcout << "|" << std::endl;


     ///*** Rubin inference
     Rcpp::List cov1L = rubin_rules_univ( cov1M, cov1_varM );
     Rcpp::List cor1L = rubin_rules_univ( cor1M, cor1_varM );

     // convert output into a matrix (set of matrices)
     Rcpp::NumericMatrix matr(VV,VV*GG);
     Rcpp::NumericMatrix matr1(VV,VV*GG);

     //---- estimated correlation
     Rcpp::NumericVector vec_pars = cor1L["pars"];
     for (int zz=0;zz<ZZ;zz++){
         int vv1=itempair_index(zz,0);
         int vv2=itempair_index(zz,1);
       for (int gg=0;gg<GG;gg++){
            matr( vv1, vv2 + gg*VV ) = vec_pars[zz*GG + gg ];
            matr( vv2, vv1 + gg*VV ) = vec_pars[zz*GG + gg];
            }
     }
     Rcpp::NumericMatrix cor1_matrix = matr;


     //---- estimated correlation
     vec_pars = cov1L["pars"];
     for (int zz=0;zz<ZZ;zz++){
         int vv1=itempair_index(zz,0);
         int vv2=itempair_index(zz,1);
       for (int gg=0;gg<GG;gg++){
            matr1( vv1, vv2 + gg*VV ) = vec_pars[zz*GG + gg ];
            matr1( vv2, vv1 + gg*VV ) = vec_pars[zz*GG + gg];
            }
     }
     Rcpp::NumericMatrix cov1_matrix = matr1;

     //*************************************************
     // OUTPUT
     return Rcpp::List::create(
         Rcpp::_["itempair_index"] = itempair_index,
         Rcpp::_["sumwgt1M"] = sumwgt1M,
         Rcpp::_["ncases1M"] = ncases1M,
         Rcpp::_["cov1"] = cov1L,
         Rcpp::_["cov1M"] = cov1M,
         Rcpp::_["cov1repM"] = cov1repM,
         Rcpp::_["cov1_varM"] = cov1_varM,
         Rcpp::_["cor1"] = cor1L,
         Rcpp::_["cor1M"] = cor1M,
         Rcpp::_["cor1repM"] = cor1repM,
         Rcpp::_["cor1_varM"] = cor1_varM,
         Rcpp::_["cor1_matrix"] = cor1_matrix,
         Rcpp::_["cov1_matrix"] = cov1_matrix
         );
}
//*************************************************************************





//*************************************************************************
//  bifie_comp_vcov_within
// [[Rcpp::export]]
Rcpp::List bifie_comp_vcov_within( Rcpp::NumericMatrix parsM,
      Rcpp::NumericMatrix parsrepM, Rcpp::NumericVector fayfac, int RR,
    int Nimp ){

    int VV = parsM.nrow();
    int NF = fayfac.size();
    double f1=0;

    //*** calculate covariance matrix for imputations
    arma::mat var_w = arma::zeros(VV,VV);
    arma::mat var_wf = arma::zeros(VV,VV*Nimp);
    Rcpp::NumericVector u_diag(VV*Nimp);

    for (int ii=0;ii<Nimp;ii++){

    for (int vv1=0;vv1<VV;vv1++){
       for (int vv2=0;vv2<VV;vv2++){
                  var_w(vv1,vv2) = 0;
       }
    }

    for (int vv1=0;vv1<VV;vv1++){
    for (int vv2=vv1;vv2<VV;vv2++){
       //@@fayfac
       f1 = fayfac[0];
       //--
    for (int rr=0;rr<RR;rr++){
          //@@fayfac
          if (NF>1){
             f1 = fayfac[rr];
                 }
          //--
       var_w(vv1,vv2) += f1 * ( parsrepM( vv1, rr+ii*RR ) - parsM( vv1, ii ) )*
            ( parsrepM( vv2, rr+ii*RR ) - parsM( vv2, ii ) );
                }
//    var_w(vv1,vv2) = fayfac[0] * var_w(vv1,vv2);
    var_w(vv2,vv1) = var_w(vv1,vv2);
                } // end vv2
    u_diag[vv1+ii*VV] = var_w(vv1, vv1);
            }  // end vv1


    for (int vv1=0;vv1<VV;vv1++){
    for (int vv2=0;vv2<VV;vv2++){
                  var_wf(vv1,vv2+ii*VV) = var_w(vv1,vv2);
                          }
                      }


    }  // end imputed dataset ii
    return Rcpp::List::create(
            Rcpp::_["u"] = var_wf,
            Rcpp::_["u_diag"] = u_diag
        );
}
//*************************************************************************



//*************************************************************************
//  bifie_comp_vcov
// [[Rcpp::export]]
Rcpp::List bifie_comp_vcov( Rcpp::NumericMatrix parsM, Rcpp::NumericMatrix parsrepM,
    Rcpp::NumericMatrix Cdes, Rcpp::NumericVector rdes, Rcpp::NumericVector Ccols,
    Rcpp::NumericVector fayfac )
{
     // number of involved variables in the test
     int VV = Ccols.size();
     int Nimp = parsM.ncol();
     double Nimp2 = Nimp + 1e-10;
     int RR = parsrepM.ncol() / Nimp;
     int df = Cdes.nrow();

     Rcpp::NumericMatrix chi2M(Nimp,2);
     arma::mat var_w = arma::zeros(VV,VV);
     arma::mat var_b = arma::zeros(VV,VV);
     Rcpp::NumericVector parsM_sel(VV);

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
         Ardes(dd,0) =  rdes[dd];
                     }

     double tmp1=0;
     double tmp2=0;

     for ( int ii=0; ii < Nimp; ii++){
         Rcpp::List res1 = bifiesurvey_rcpp_wald_test_vcov(  VV,  Ccols, parsM, parsrepM,
              ii,  RR, fayfac, ACdes,  Ardes );

         Rcpp::NumericMatrix chi2a=res1["chi2"];
         chi2M(ii,0) = chi2a(0,0);
         tmp1 += chi2M(ii,0);
         chi2M(ii,1) = sqrt( chi2M(ii,0) );
         tmp2 += chi2M(ii,1);
         Rcpp::NumericMatrix var_w1 = res1["var_w"];
         for (int vv1=0;vv1<VV;vv1++){
            for (int vv2=0;vv2<VV;vv2++){
                  var_w(vv1,vv2) += var_w1(vv1,vv2);
                              }
                          }
             }

     // calculate ARIV
     double eps=1e-10;
     double ariv = tmp1 - Nimp * pow( tmp2 / Nimp, 2.0 );
     ariv = ariv / ( Nimp - 1 + eps ) * ( 1 + 1 / Nimp2 );

     // calculate D2 statistic
     double D2 = tmp1 / Nimp2 - (Nimp2+1)/(Nimp2-1+eps) * ariv;
     D2 = D2 / ( 1 + ariv );
     // calculate degrees of freedom
     double df2 = df;
     double nu3 = 1000;
     if ( Nimp > 1 ){
       nu3 = pow( df2, - 3 / Nimp2 ) * ( Nimp2 - 1 ) *
                 pow( 1 + 1 / ( ariv + eps ), 2 );
                   }

     // double p_D2 = Rf_pf( D2, df, nu3, FALSE, FALSE );

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
     // parsM( Ccols[vv1], ii ) )

     for (int vv1=0;vv1<VV;vv1++){
     for (int vv2=0;vv2<VV;vv2++){
       for ( int ii=0; ii<Nimp; ii++){
           var_b(vv1,vv2) += ( parsM( Ccols[vv1], ii ) - parsM_sel[vv1] ) *
                  ( parsM( Ccols[vv2], ii ) - parsM_sel[vv2] );
                     }
           var_b(vv1,vv2) = 1 / ( Nimp2 - 1 ) * var_b( vv1, vv2 );
                   }
               }

     arma::mat ariv_D1a = arma::mat( var_b  );
     double ariv_D1 = 0;
     for (int vv=0;vv<VV;vv++){
         ariv_D1 += ariv_D1a(vv,vv);
                 }
     ariv_D1 = ariv_D1 * ( 1 + 1 / Nimp2 ) / df;
     arma::mat var_t1 = arma::mat( (1+ariv_D1) * var_w );

     // hypothesis matrix
     arma::mat var_hyp = arma::mat( ACdes * var_t1 * arma::trans( ACdes) );
     // compute inverse of variance matrix of hypothesis
     arma::mat var_hypinv = arma::mat( var_hyp );
     // parameter vector
     arma::colvec parm_vec= arma::zeros(VV,1);
     for (int vv=0;vv<VV;vv++){
             parm_vec(vv,0) = parsM_sel[vv];
                     }
     // hypothesis statistic
     arma::mat hyp_stat = arma::mat( ACdes * parm_vec - Ardes );
     arma::mat D1 = arma::mat( arma::trans( hyp_stat ) * var_hypinv * hyp_stat );
     // D1(0,0) = D1(0,0) / df;
     // according to Enders (2010, p. 236), D1 must be divided by df,
     // but this is (could be) an error?

     // calculate nu2
     double nu2 = 1 + ( 1 - 2 / ( df * Nimp2 - df ) * 1 / ariv_D1 );
     nu2 = 4 + ( df * Nimp2 - df - 4 ) * nu2 * nu2;

     // double tmp11 = D1(0,0);

     // double p_D1 = Rf_pf( tmp11, df, nu2, FALSE, FALSE );

     //*************************************************
     // OUTPUT

     return Rcpp::List::create(
         Rcpp::_["chi2M"] = chi2M,
         Rcpp::_["ariv"] = ariv,
         Rcpp::_["D2"] = D2,
         Rcpp::_["df"] = df,
         Rcpp::_["nu2"] = nu2,
         Rcpp::_["nu3"]=nu3,
         Rcpp::_["Nimp"] = Nimp,
         Rcpp::_["RR"] = RR,
         Rcpp::_["fayfac"] = fayfac,
         Rcpp::_["var_w"] = var_w,
         Rcpp::_["var_b"] = var_b,
         Rcpp::_["D1"] = D1,
         Rcpp::_["Ccols"] = Ccols,
         Rcpp::_["parsM_sel"] = parsM_sel
         );
}

//*************************************************************************

//*************************************************************************
//  bifie_test_univar
// [[Rcpp::export]]
Rcpp::List bifie_test_univar( Rcpp::NumericMatrix mean1M, Rcpp::NumericMatrix sd1M,
          Rcpp::NumericMatrix sumweightM, int GG, Rcpp::NumericVector group_values,
          Rcpp::NumericMatrix mean1repM, Rcpp::NumericMatrix sd1repM,
          Rcpp::NumericMatrix sumweightrepM, Rcpp::NumericVector fayfac  ){

     int HH = mean1M.nrow();
     int VV = HH / GG;
     int GG2 = GG * (GG-1) / 2;
     int Nimp = sd1M.ncol();
     int RR = sd1repM.ncol() / Nimp;
     Rcpp::NumericMatrix dstatM(VV*GG2,Nimp);
     Rcpp::NumericMatrix dstatrepM1(VV*GG2,Nimp*RR);
     Rcpp::NumericMatrix dstat_varM(VV*GG2,Nimp);
     Rcpp::NumericMatrix eta2M(VV,Nimp);
     Rcpp::NumericMatrix eta2repM1(VV,Nimp*RR);
     Rcpp::NumericMatrix eta2_varM(VV,Nimp);
     Rcpp::NumericVector eta2V(1);

     // matrix of group values
     Rcpp::NumericMatrix group_values_matrix(GG2,2);
     int ii=0;
     for ( int gg1=0; gg1 < GG - 1; gg1++){
     for (int gg2=gg1+1; gg2 < GG; gg2++){
         group_values_matrix(ii,0) = group_values[gg1];
         group_values_matrix(ii,1) = group_values[gg2];
                ii++;
                 }
             }
     Rcpp::Rcout << "|";

     // loop over imputations
     for ( int ii=0; ii < Nimp; ii++){
     for (int vv=0; vv < VV; vv++){

     // dataset ii
     Rcpp::NumericMatrix mean1M_ii = mean1M( Rcpp::Range(vv*GG,vv*GG+GG-1), Rcpp::Range(ii,ii) );
     Rcpp::NumericMatrix sd1M_ii = sd1M( Rcpp::Range(vv*GG,vv*GG+GG-1), Rcpp::Range(ii,ii) );
     Rcpp::NumericMatrix sumweightM_ii = sumweightM( Rcpp::Range(0,GG-1), Rcpp::Range(ii,ii) );

     Rcpp::List res = bifiehelpers_etasquared( mean1M_ii, sd1M_ii, sumweightM_ii,  GG );
     eta2V = matr2vec( res["eta2"] );
     Rcpp::NumericVector dstatV = matr2vec( res["dstat"] );

     // analysis replicate weights
     Rcpp::NumericMatrix mean1M_rr = mean1repM( Rcpp::Range(vv*GG,vv*GG+GG-1), Rcpp::Range(ii*RR,ii*RR + RR-1) );
     Rcpp::NumericMatrix sd1M_rr = sd1repM( Rcpp::Range(vv*GG,vv*GG+GG-1), Rcpp::Range(ii*RR,ii*RR + RR-1)  );
     Rcpp::NumericMatrix sumweightM_rr = sumweightrepM( Rcpp::Range(0,GG-1), Rcpp::Range(ii*RR,ii*RR + RR-1)  );
     Rcpp::List res1 = bifiehelpers_etasquared( mean1M_rr, sd1M_rr, sumweightM_rr,  GG );
     Rcpp::NumericMatrix eta2repM = res1["eta2"];
     Rcpp::NumericMatrix dstatrepM = res1["dstat"];

     // compute standard errors
     // eta squared
     Rcpp::NumericVector eta2_var = varjack_helper( eta2V, eta2repM, fayfac );
     // d statistics
     Rcpp::NumericVector dstat_var = varjack_helper( dstatV, dstatrepM, fayfac );
     // adjusted eta squared statistic
     // Rcpp::List res3 = varjack_bias_helper( eta2V, eta2repM, fayfac );
     // Rcpp::NumericVector eta2adj_tmp = res3["pars_bias"];
     // Rcpp::NumericVector eta2adj_var = res3["pars_var"];
     // double eta2adj = eta2V[0] + (RR-1)*( eta2adj_tmp[0] - eta2V[0] );
     // double eta2adj = RR * eta2V[0] - (RR-1)* eta2adj_tmp[0];

     eta2M(vv,ii) = eta2V[0];
     eta2_varM(vv,ii) = eta2_var[0];
     for (int zz=0;zz<GG2;zz++){
         dstatM(zz+vv*GG2,ii) = dstatV[zz];
         dstat_varM(zz+vv*GG2,ii) = dstat_var[zz];
         for ( int rr=0;rr<RR;rr++){
             dstatrepM1(zz+vv*GG2,rr+ii*RR) = dstatrepM(zz,rr);
                     }
                 }

     //    int GG2 = GG * (GG - 1 ) / 2;
     //    Rcpp::NumericMatrix dstat(GG2, WW );
     for ( int rr=0;rr<RR;rr++){
              eta2repM1(vv,rr+ii*RR) = eta2repM(0,rr);
                          }


             } // end vv
     Rcpp::Rcout << "-" <<  std::flush;

          }  // end ii;  end multiple imputations

     Rcpp::Rcout << "|" << std::endl;
     //----------------------

     ///*** Rubin inference
     Rcpp::List eta2L = rubin_rules_univ( eta2M, eta2_varM );
     Rcpp::List dstatL = rubin_rules_univ( dstatM, dstat_varM );

     //*************************************************
     // OUTPUT
     return Rcpp::List::create(
         Rcpp::_["eta2L"] = eta2L,
         Rcpp::_["eta2M"] = eta2M,
         Rcpp::_["eta2repM"] = eta2repM1,
         Rcpp::_["eta2_varM"] = eta2_varM,
         Rcpp::_["dstatL"] = dstatL,
         Rcpp::_["dstatM"] = dstatM,
         Rcpp::_["dstatrepM"] = dstatrepM1,
         Rcpp::_["dstat_varM"] = dstat_varM,
         Rcpp::_["group_values_matrix"] = group_values_matrix
         );
}
//*************************************************************************

//*************************************************************************
//  bifie_crosstab
// [[Rcpp::export]]
Rcpp::List bifie_crosstab( Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1,
     Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector vars_values1,
     Rcpp::NumericVector vars_index1, Rcpp::NumericVector vars_values2,
     Rcpp::NumericVector vars_index2, Rcpp::NumericVector fayfac,
     Rcpp::NumericVector NI, Rcpp::NumericVector group_index1,
     Rcpp::NumericVector group_values ){

     int Nimp = NI[0];
     int RR = wgtrep.ncol();
     int N = wgt1.nrow();
     int VV1 = vars_values1.size();
     int VV2 = vars_values2.size();
     int NV = datalist.ncol();
     int GG=group_values.size();

     Rcpp::NumericMatrix dat1(N,NV);
     int ZZ = VV1*VV2*GG;

     // design matrix
     Rcpp::NumericMatrix design_pars(ZZ,5);
     int zz=0;
     for (int gg=0;gg<GG;gg++){
     for (int vv1=0;vv1<VV1;vv1++){
     for (int vv2=0;vv2 <VV2;vv2++){
         design_pars(zz,0) = vars_values1[vv1];
         design_pars(zz,1) = vars_values2[vv2];
         design_pars(zz,2) = group_values[gg];
         design_pars(zz,3) = vv1;
         design_pars(zz,4) = vv2;
         zz++;
         }  // end vv2
     }  // end vv1
     }  // end gg

     int CTP = 3*ZZ + VV1*GG + VV2*GG + 2*GG + GG + 3*GG + 3*GG;
     Rcpp::NumericMatrix ncasesM( ZZ, Nimp );
     Rcpp::NumericMatrix ncases_ggM( GG, Nimp );
     Rcpp::NumericMatrix sumwgtM( ZZ, Nimp);
     Rcpp::NumericMatrix ctparsM( CTP, Nimp);
     Rcpp::NumericMatrix ctpars_varM( CTP, Nimp);
     Rcpp::NumericMatrix ctparsrepM( CTP, Nimp*RR);

     Rcpp::Rcout << "|";

     ///****** loop imputed datasets

     for (int ii=0; ii <Nimp; ii++){

         dat1 = datalist( Rcpp::Range( ii*N+0, ii*N+ (N-1) ), Rcpp::Range(0,NV-1) );

         //*** analysis for original data
         Rcpp::List res20 = bifiehelpers_crosstab( dat1,  wgt1, group_values, group_index1,
             vars_values1,  vars_index1, vars_values2, vars_index2, design_pars );
         Rcpp::NumericMatrix ncases = res20["ncases"];
         Rcpp::NumericMatrix ncases_gg = res20["ncases_gg"];
         Rcpp::NumericMatrix sumwgt = res20[ "sumwgt"];
         Rcpp::NumericMatrix sumwgt_gg = res20["sumwgt_gg"];
         Rcpp::NumericVector ctpars = matr2vec(res20["crosstab_pars"]);

         //*** analysis for replicated datasets
         Rcpp::List res3 = bifiehelpers_crosstab( dat1,  wgtrep, group_values, group_index1,
             vars_values1,  vars_index1, vars_values2, vars_index2, design_pars );
         Rcpp::NumericMatrix ctparsrep = res3["crosstab_pars"];

         // compute standard errors
         Rcpp::NumericVector ctpars_var = varjack_helper( ctpars, ctparsrep, fayfac );

         // data handling
         ncasesM(_,ii) = ncases(_,0);
         ncases_ggM(_,ii) = ncases_gg(_,0);
         sumwgtM(_,ii) = sumwgt(_,0);
         ctparsM(_,ii) = ctpars;
         ctpars_varM(_,ii) = ctpars_var;
         for (int rr=0;rr<RR;rr++){
             ctparsrepM(_,rr+ii*RR) = ctparsrep(_,rr);
                     }

         Rcpp::Rcout << "-" <<  std::flush;
             }

     Rcpp::Rcout << "|" << std::endl;

     ///*** Rubin inference
     Rcpp::List ctparsL = rubin_rules_univ( ctparsM, ctpars_varM );

     //*************************************************
     // OUTPUT
     return Rcpp::List::create(
         Rcpp::_["design_pars"] = design_pars,
         Rcpp::_["ncases_ggM"] = ncases_ggM,
         Rcpp::_["ncasesM"] = ncasesM,
         Rcpp::_["sumwgtM"] = sumwgtM,
         Rcpp::_["ctparsL"] = ctparsL,
         Rcpp::_["ctparsM"] = ctparsM,
         Rcpp::_["ctparsrepM"] = ctparsrepM,
         Rcpp::_["ctpars_varM"] = ctpars_varM
         );
}
//*************************************************************************

//*************************************************************************
//  bifie_by
// [[Rcpp::export]]
Rcpp::List bifie_by( Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1,
         Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector vars_index,
         Rcpp::NumericVector fayfac, Rcpp::NumericVector NI,
         Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values,
         Rcpp::Function userfct ){

     int Nimp = NI[0];
     int N = wgt1.nrow();
     int VV = vars_index.size();
     int NV = datalist.ncol();
     int group_index = group_index1[0];
     int GG=group_values.size();
     Rcpp::NumericMatrix dat1(N,NV);
     int WW = wgtrep.ncol();

     // start with a calculation to compute the number of parameters
     Rcpp::NumericVector w = wgt1(_,0);
     int ii=0;
     dat1 = datalist( Rcpp::Range( ii*N+0, ii*N+ (N-1) ), Rcpp::Range(0,NV-1) );
     Rcpp::NumericMatrix X(N,VV);
     for (int vv=0;vv<VV;vv++){
         X(_,vv) = dat1(_, vars_index[vv] );
     }
     Rcpp::NumericVector pars= userfct(X,w);
     int NP = pars.size();


     Rcpp::NumericMatrix ncasesM(GG,Nimp);
     Rcpp::NumericMatrix sumwgtM(GG,Nimp);

     Rcpp::NumericMatrix parsM(NP*GG,Nimp);
     Rcpp::NumericMatrix pars_varM(NP*GG,Nimp);
     Rcpp::NumericVector pars1(NP*GG);
     Rcpp::NumericMatrix pars1rep(NP*GG,WW);
     Rcpp::NumericMatrix pars1repM(NP*GG,WW*Nimp);
     Rcpp::NumericVector pars1_var(NP*GG);

     Rcpp::Rcout << "|";

     // dataset ii
     for (ii=0;ii<Nimp;ii++){  // beg ii
         dat1 = datalist( Rcpp::Range( ii*N+0, ii*N+ (N-1) ), Rcpp::Range(0,NV-1) );

         //-- compute dimensions
         for (int nn=0;nn<N;nn++){ // beg nn
         for (int gg =0;gg<GG;gg++){  // beg gg
         if ( dat1(nn,group_index) == group_values[gg] ){
             ncasesM(gg,ii) ++;
             sumwgtM(gg,ii) += wgt1(nn,0);
             break;
                     }
                 }  // end gg
             } // end nn

         int hh=0;

         for ( int gg=0; gg <GG; gg++){  // beg gg
            //-- evaluate function    for original dataset ii
             Rcpp::NumericMatrix X1( ncasesM(gg,ii),VV);
             Rcpp::NumericVector w1( ncasesM(gg,ii) );
             hh=0;
             for (int nn=0;nn<N;nn++){  // beg nn
             if ( dat1(nn, group_index ) == group_values[gg] ){ // beg group val
                 for (int vv=0;vv<VV;vv++){
                     X1(hh,vv) = dat1( nn, vars_index[vv] );
                     }
                 w1[hh] = wgt1(nn,0);
                 hh++;
                 }  // end if group val
             }  // end nn

             Rcpp::NumericVector pars_res = userfct( X1, w1 );
             for (int pp=0;pp<NP;pp++){  // beg pp
                 parsM( pp + gg*NP, ii ) = pars_res[pp];
                 pars1[ pp + gg*NP ] = pars_res[pp];
                         }  // end pp
             //*** evaluate user function for replicated datasets
             for (int ww=0;ww<WW;ww++){  // beg ww
                 Rcpp::NumericVector w2( ncasesM(gg,ii) );
                 hh=0;
                 for (int nn=0;nn<N;nn++){  // beg nn
                  if ( dat1(nn, group_index ) == group_values[gg] ){
                     w2[hh] = wgtrep(nn,ww);
                     hh++;
                     }
                 }  // end nn
                 Rcpp::NumericVector pars_res2 = userfct( X1, w2 );
                 for (int pp=0;pp<NP;pp++){  // beg pp
                     pars1repM( pp + gg*NP, ww + ii*WW ) = pars_res2[pp];
                     pars1rep( pp + gg*NP, ww ) = pars_res2[pp];
                             }  // end pp
                     } // end ww
             }  // end gg

         // compute standard errors
         pars1_var = varjack_helper( pars1, pars1rep, fayfac );
         pars_varM(_,ii) = pars1_var;
     Rcpp::Rcout << "-" <<  std::flush;

          }  // end ii;  end multiple imputations

     Rcpp::Rcout << "|" << std::endl;

     ///*** Rubin inference
     Rcpp::List parsL = rubin_rules_univ( parsM, pars_varM );
     //*************************************************
     // OUTPUT
     return Rcpp::List::create(
         Rcpp::_["WW"] = WW,
         Rcpp::_["N"] = N,
         Rcpp::_["NP"] = NP,
         Rcpp::_["userfct"] = userfct,
         Rcpp::_["sumwgtM"] = sumwgtM,
         Rcpp::_["parsrepM"] = pars1repM,
         Rcpp::_["parsM"] = parsM,
         Rcpp::_["pars_varM"] = pars_varM,
         Rcpp::_["ncasesM"] = ncasesM,
         Rcpp::_["parsL"] = parsL
         );
}
//*************************************************************************

//*************************************************************************
//  bifie_hist
// [[Rcpp::export]]
Rcpp::List bifie_hist( Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1,
          Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector vars_index,
          Rcpp::NumericVector fayfac, Rcpp::NumericVector NI,
          Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values,
          Rcpp::NumericVector breaks ){

     int Nimp = NI[0];
     int N = wgt1.nrow();
     int NV = datalist.ncol();
     int group_index = group_index1[0];
     int GG=group_values.size();
     int BB=breaks.size() - 1;
     Rcpp::NumericMatrix dat1(N,NV);
     Rcpp::NumericMatrix countsM(BB*GG,Nimp);
     Rcpp::NumericMatrix sumwgtM(BB*GG,Nimp);
     Rcpp::NumericMatrix ncasesM(GG,Nimp);
     Rcpp::NumericMatrix sumwgt_ggM(GG,Nimp);
     Rcpp::NumericVector counts(BB*GG);
     Rcpp::NumericVector sumwgt(BB*GG);
     Rcpp::NumericVector sumwgt_gg(GG);

     int bb=0;

     Rcpp::Rcout << "|";

     for ( int ii=0; ii < Nimp; ii++ ){

     dat1 = datalist( Rcpp::Range( ii*N+0, ii*N+ (N-1) ), Rcpp::Range(0,NV-1) );

     for (int nn=0; nn<N;nn++){ // beg nn
     for (int gg=0; gg<GG;gg++){ // beg gg
     if ( dat1(nn, group_index ) == group_values[gg] ){  // beg if group val
        if ( ! R_IsNA( dat1( nn, vars_index[0] ) ) ){  // beg non NA
         ncasesM(gg,ii) ++;
         sumwgt_ggM(gg,ii) += wgt1(nn,0);
         bb = 0;
             while (bb < BB ){  // beg while bb
                  if ( dat1( nn, vars_index[0] ) < breaks[bb+1] ){ // beg if larger
                          countsM( bb + gg*BB, ii ) ++;
                          sumwgtM( bb + gg*BB, ii ) += wgt1(nn,0);
                          bb = BB;
                                  } // end if larger
                   bb ++;
                     } // end while bb
             break;
            } // if non NA
     }  // end if group val
     } // end if gg
     } // end nn

     for (int hh=0; hh < BB*GG; hh++){
         counts[hh] += countsM(hh,ii);
         sumwgt[hh] += sumwgtM(hh,ii);
                     }
     for (int hh=0; hh < GG; hh++){
         sumwgt_gg[hh] += sumwgt_ggM(hh,ii);
                     }

     Rcpp::Rcout << "-" <<  std::flush;

     } // end ii

     Rcpp::Rcout << "|" << std::endl;


     //*********************************
     // compute statistics
     for (int hh=0; hh < BB*GG; hh++){
         counts[hh] = counts[hh] / Nimp;
         sumwgt[hh] = sumwgt[hh] / Nimp;
     }
     for (int hh=0; hh < GG; hh++){
         sumwgt_gg[hh] = sumwgt_gg[hh] / Nimp;
     }

     // mid points
     Rcpp::NumericVector mids(BB-1);
     for (int bb=0;bb<BB-1;bb++){
         mids[bb] = ( breaks[bb] + breaks[bb+1] ) / 2.0;
     }

     // density
     Rcpp::NumericVector density_vec(BB*GG);
     Rcpp::NumericVector relfreq(BB*GG);

     for (int gg=0;gg<GG; gg++){
        for (int bb=0;bb<BB;bb++){
           relfreq[bb+gg*BB] = sumwgt[bb+gg*BB] / sumwgt_gg[gg];
           density_vec[bb+gg*BB] = relfreq[bb+gg*BB] / ( breaks[bb+1] - breaks[bb] );
        }
     }
     //*************************************************
     // OUTPUT
     return Rcpp::List::create(
          Rcpp::_["BB"] = BB,
         Rcpp::_["breaks"] = breaks,
         Rcpp::_["mids"] = mids,
         Rcpp::_["sumwgtM"] = sumwgtM,
         Rcpp::_["countsM"] = countsM,
         Rcpp::_["ncasesM"] = ncasesM,
         Rcpp::_["counts"] = counts,
         Rcpp::_["sumwgt"] = sumwgt,
         Rcpp::_["sumwgt_ggM"] = sumwgt_ggM,
         Rcpp::_["sumwgt_gg"] = sumwgt_gg,
         Rcpp::_["relfreq"] = relfreq,
         Rcpp::_["density_vec"] = density_vec
         );
}
//*************************************************************************

//*************************************************************************
//  bifie_ecdf
// [[Rcpp::export]]
Rcpp::List bifie_ecdf( Rcpp::NumericMatrix datalist, Rcpp::NumericMatrix wgt1,
       Rcpp::NumericMatrix wgtrep, Rcpp::NumericVector vars_index,
       Rcpp::NumericVector fayfac, Rcpp::NumericVector NI,
       Rcpp::NumericVector group_index1, Rcpp::NumericVector group_values,
       Rcpp::NumericVector breaks, int quanttype, int maxval ){

     int Nimp = NI[0];
     int N = wgt1.nrow();
     int VV = vars_index.size();
     int NV = datalist.ncol();
     int GG=group_values.size();
     int BB=breaks.size();
     int ZZ=VV*GG*BB;
     Rcpp::NumericMatrix ecdfM(ZZ,Nimp);
     Rcpp::NumericVector ecdfMtemp(ZZ);
     Rcpp::NumericMatrix ncasesM(VV*GG,Nimp);
     Rcpp::NumericMatrix sumwgtM(VV*GG,Nimp);
     Rcpp::NumericMatrix dat1(N,NV);

     for (int ii=0;ii<Nimp;ii++){ // beg dataset ii

     // int ii=0;
         dat1 = datalist( Rcpp::Range( ii*N+0, ii*N+ (N-1) ), Rcpp::Range(0,NV-1) );

         ecdfMtemp = bifie_helper_ecdf( dat1, wgt1, breaks,
              group_values, group_index1,vars_index, ii,
              ncasesM,  sumwgtM, maxval, quanttype );

         for (int zz=0;zz<ZZ;zz++){
             ecdfM(zz,ii) = ecdfMtemp[zz];
                     }
     // Rcpp::Rcout << "ii= " <<  ii <<  std::flush << std::endl;
         } // end dataset ii

     //********************************
     // average ecdf
     Rcpp::NumericVector ecdf(ZZ);
     for (int zz=0;zz<ZZ;zz++){
       for (int ii=0;ii<Nimp;ii++){
         ecdf[zz] += ecdfM(zz,ii);
                 }
         ecdf[zz] = ecdf[zz] / Nimp;
         }

     //*************************************************
     // OUTPUT
     return Rcpp::List::create(
         Rcpp::_["BB"] = BB,
         Rcpp::_["yval"] = breaks,
         Rcpp::_["sumwgtM"] = sumwgtM,
         Rcpp::_["ncasesM"] = ncasesM,
         Rcpp::_["ecdfM"] = ecdfM,
         Rcpp::_["ecdf"] = ecdf
         );
}
//*************************************************************************


//*************************************************************************
//  bifie_fasttable
// [[Rcpp::export]]
Rcpp::List bifie_fasttable( Rcpp::NumericMatrix datavec ){

     int N = datavec.nrow();
     arma::colvec vals_temp(N);
     int ii=0;
     for (int nn=0;nn<N;nn++){
         if ( ! R_IsNA( datavec(nn,0) ) ){
             vals_temp(ii,0) = datavec(nn,0);
             ii ++;
                         }
                     }

     int N1 = ii-1;

     // arma::mat vec_unique = unique(vals_temp);
     arma::mat vec_sort = arma::sort( vals_temp( arma::span(0,N1), arma::span(0,0) ) );

     // create result table
     Rcpp::NumericMatrix tableM(N1,2);
     ii = 0;
     tableM(ii,0) = vec_sort(0,0);
     tableM(ii,1) = 1;

     for (int nn=1;nn<N1+1;nn++){
         if (vec_sort(nn,0) == tableM(ii,0) ){
             tableM(ii,1) ++;
                 } else {
             ii ++;
                 tableM(ii,0) = vec_sort(nn,0);
                 tableM(ii,1) = 1;
                     }
     }

     //*************************************************
     // OUTPUT
     return Rcpp::List::create(
         Rcpp::_["vec_sort"] = vec_sort,
         Rcpp::_["tableM"] = tableM,
         Rcpp::_["N_unique"] = ii+1
         );
}

//*************************************************************************

//*************************************************************************
//  bifie_table1_character
// [[Rcpp::export]]
Rcpp::List bifie_table1_character( Rcpp::CharacterVector datavec ){

     int N = datavec.size();
     Rcpp::CharacterVector uii = Rcpp::unique( datavec );
     Rcpp::IntegerVector indii = Rcpp::match( datavec, uii );
     int Nval = uii.size();
     Rcpp::NumericVector tableM(Nval);

     for (int nn=0;nn<N;nn++){
         tableM[ indii[nn] - 1 ] ++;
     }
     //*************************************************
     // OUTPUT
     return Rcpp::List::create(
         Rcpp::_["table_names"] = uii,
         Rcpp::_["tableM"] = tableM
         );
}

//*************************************************************************

//*************************************************************************
//  bifie_mla2
// [[Rcpp::export]]
Rcpp::List bifie_mla2( Rcpp::NumericMatrix X_list, Rcpp::NumericMatrix Z_list,
       Rcpp::NumericVector y_list, Rcpp::NumericVector wgttot,
       Rcpp::NumericVector wgtlev2, Rcpp::NumericVector wgtlev1, double globconv,
       int maxiter, Rcpp::NumericVector group, Rcpp::NumericVector group_values,
       Rcpp::NumericVector cluster, Rcpp::NumericMatrix wgtrep,
       int Nimp, Rcpp::NumericVector fayfac,
       Rcpp::NumericMatrix recov_constraint, int is_rcov_constraint ){

     // new declarations
     int NZ = Z_list.ncol();
     int NX = X_list.ncol();
     int N = wgttot.size();
     int NC = wgtlev2.size();
     int GG = group_values.size();
     int RR = wgtrep.ncol();
     //*** number of constraints
     int NRC = recov_constraint.nrow();

     double eps = 1E-8;

     Rcpp::NumericMatrix idcluster_table2;
     Rcpp::NumericVector pars;

     // estimated parameters
     int NP = NX + NZ*NZ + 1 + 13;
     int NPtot = NP * GG;
     Rcpp::NumericMatrix parsM(NPtot, Nimp );
     Rcpp::NumericMatrix parsVar(NPtot, Nimp );
     Rcpp::NumericMatrix parsMrep(NPtot, RR*Nimp );
     Rcpp::NumericMatrix iterM( GG, Nimp );
     Rcpp::NumericMatrix iterMrep( GG, RR*Nimp );
     Rcpp::NumericMatrix fvcovM( GG*NX, NX*Nimp );
     Rcpp::NumericVector Npers(GG);
     Rcpp::NumericVector Nclusters(GG);
     Rcpp::NumericVector pars0(NP);
     Rcpp::NumericMatrix pars20(NP,RR);

     Rcpp::NumericMatrix Sigma_W_yXM(Nimp*NX,NX);
     Rcpp::NumericMatrix Sigma_B_yXM(Nimp*NX,NX);
     Rcpp::NumericMatrix Sigma_W_yZM(Nimp*NZ,NZ);
     Rcpp::NumericMatrix Sigma_B_yZM(Nimp*NZ,NZ);
     Rcpp::NumericMatrix totmean_yXM(Nimp,NX);
     Rcpp::NumericMatrix totmean_yZM(Nimp,NZ);

     // design matrices
     Rcpp::NumericVector y(N);
     Rcpp::NumericMatrix Z(N,NZ);
     Rcpp::NumericMatrix X(N,NX);


     int maxiter_rep = maxiter;
     if (RR==1){ maxiter_rep = 1; }

     /////////////////////////////////////////////
     // imputations
     /////////////////////////////////////////////

     for (int imp=0; imp<Nimp; imp++){
     // int imp=0; // imputation imp

     int vv=0;

     for (int ii=0;ii<N;ii++){
         y[ii] = y_list[ ii + imp*N ];
         for (int jj=0;jj<NX;jj++){
                 X( ii, jj ) = X_list( ii + imp* N, jj );
                            }
         for (int jj=0;jj<NZ;jj++){
                 Z( ii, jj ) = Z_list( ii + imp* N, jj );
                            }
                 }

     //************************************
     // group-wise analysis and missings (first imputed dataset)
     Rcpp::List res41 = create_dummies_mla2( GG, group, X, Z, y );
     Rcpp::NumericMatrix dummy_inds = res41["dummy_inds"];
     Rcpp::NumericVector N_group = res41["N_group"];
     if (imp==0){
     for (int gg=0;gg<GG;gg++){
         Npers[gg] = N_group[gg];
                 }
             }

     // collect data for group gg
     // X, Z, y, wgtlev1, wgttot, wgtlev2, cluster, group
     //Rcpp::NumericMatrix pars2;
     Rcpp::NumericVector iter3;
     Rcpp::List res36;

     //--- loop group
     for (int gg=0;gg<GG;gg++){

     Rcpp::NumericVector pars;
     //Rcpp::NumericMatrix fvcov;
     int iter2=0;

     int N__ = N_group[gg];
     Rcpp::NumericMatrix X__(N__,NX);
     Rcpp::NumericMatrix Z__(N__,NZ);
     Rcpp::NumericMatrix wgtrep__(N__,RR);
     Rcpp::NumericVector y__(N__);
     Rcpp::NumericVector wgtlev1__(N__);
     Rcpp::NumericVector wgttot__(N__);
     Rcpp::NumericVector cluster__(N__);
     Rcpp::NumericVector group__(N__);
     Rcpp::NumericVector clustertable_temp(NC);

     int kk=0;

     for( int nn=0; nn < N; nn++){
         if ( dummy_inds(nn,gg) == 1){
                for (int ii=0;ii<NX;ii++){
                X__(kk,ii) = X(nn,ii);
                     }
                for (int ii=0;ii<NZ;ii++){
                Z__(kk,ii) = Z(nn,ii);
                     }
                for (int ii=0;ii<RR;ii++){
                wgtrep__(kk,ii) = wgtrep(nn,ii);
                     }
             y__[kk] = y[nn];
             wgtlev1__[kk] = wgtlev1[nn];
             wgttot__[kk] = wgttot[nn];
             cluster__[kk] = cluster[nn];
             group__[kk] = group[nn];
             kk ++;
             }
         }

     //** count clusters
     int NC__ = 0;
     double cl = -999;
     for ( int nn = 0; nn < N__; nn++){
         if ( cluster__[nn] > cl ){
             cl = cluster__[nn];
             clustertable_temp[ NC__ ] = cl;
             NC__ ++;
                     }
                 }

     Rcpp::NumericVector wgtlev2__(NC__);
     cl=0;
     for (int nn=0;nn<NC__;nn++){
           wgtlev2__[nn] = wgtlev2[ clustertable_temp[nn] ];
           cl += wgtlev2__[nn];
                       }

     for (int nn=0;nn<NC__;nn++){
           wgtlev2__[nn] = wgtlev2__[nn] * NC__ / cl;
                       }
     for (int nn=0;nn<N__;nn++){
           wgttot__[nn] = wgttot__[nn] * NC__ / cl;
           if (wgttot__[nn] < eps ){
                     wgttot__[nn] = eps;
                 }
                       }

     if (imp == 0 ){
         Nclusters[gg] = NC__;
             }

     //    Rcpp::Rcout << "a150" <<    std::flush << std::endl;

     //***********************************
     // create cluster table
     idcluster_table2 = create_idclustertable( group__, cluster__,  NC__);
     //    Rcpp::Rcout << "a200" <<    std::flush << std::endl;

     //**********************************
     // get initial values (if iter=0);
     // arma::mat Xa(X.begin(), N, NX, false);
     arma::mat Xa = rcppmat2armamat( X__ )["armamat"];
     Rcpp::List res1 = mla2_inits( Xa, X__, Z__, y__, NZ, wgtlev1__, wgttot__ );
     arma::mat theta_init = res1["theta"];
     arma::mat Tmat_init = res1["Tmat"];
     arma::mat sig2_init = res1["sig2"];
     if (imp > 0 ){
       vv = 0;
       // theta inits
       for (int ii=0;ii<NX;ii++){
             theta_init(vv,0) = parsM( vv + gg*NP, 0 );
             vv ++;
          }
       // Tmat inits
       for (int ii=0;ii<NZ;ii++){
            for (int jj=0;jj<NZ;jj++){
                 Tmat_init(ii,jj) = parsM(vv + gg*NP, 0 );
                 if ( jj>ii){
                      Tmat_init(jj,ii) = Tmat_init(ii,jj);
             }
                 vv ++;
            }
          }
        // sig2 inits
        sig2_init(0,0) = parsM(vv+gg*NP,0);
    }
     //-------------- DATASET ORIGINAL --------------------------

     // print progress
     Rcpp::Rcout << " " << std::endl;
     Rcpp::Rcout << "Imputation " <<  imp+1 <<
        " | Group " << gg+1 <<   " |" << std::flush;

     //**********************************
     // rescaling weights
     Rcpp::NumericVector wgtlev1a = rescale_lev1weights( idcluster_table2,
                          wgtlev1__ );

     //*********************************
     // estimate model (FIML)
     Rcpp::List res32 = bifie_mla2_estimation( theta_init, Tmat_init,
          sig2_init, NX, NZ,  NC__,  N__, X__,  Z__, y__, wgtlev2__,
          wgtlev1a, wgttot__, idcluster_table2, globconv, maxiter,
       recov_constraint, is_rcov_constraint, NRC );
     pars = res32["pars_"]; // extract estimated parameters
     for (int pp=0;pp<NP;pp++){
        parsM( pp + gg*NP, imp ) = pars[pp];
                    }
     Rcpp::NumericMatrix fvcov = res32["fvcov"];
     for (int ii=0;ii<NX;ii++){
     for (int jj=0;jj<NX;jj++){
          fvcovM( ii + gg*NX, jj + NX*imp ) = fvcov(ii,jj);
                 }
             }
     // save covariance decompositions;
     Rcpp::List res32a = res32["postproc"];
     Rcpp::NumericMatrix Sigma_W_yX = res32a["Sigma_W_yX"];
     Rcpp::NumericMatrix Sigma_B_yX = res32a["Sigma_B_yX"];
     Rcpp::NumericMatrix Sigma_W_yZ = res32a["Sigma_W_yZ"];
     Rcpp::NumericMatrix Sigma_B_yZ = res32a["Sigma_B_yZ"];
     Rcpp::NumericVector totmean_yX = res32a["totmean_yX"];
     Rcpp::NumericVector totmean_yZ = res32a["totmean_yZ"];

     for (int ii=0;ii<NX;ii++){
     for (int jj=0;jj<NX;jj++){
         Sigma_W_yXM(ii+NX*imp,jj) = Sigma_W_yX(ii,jj);
         Sigma_B_yXM(ii+NX*imp,jj) = Sigma_B_yX(ii,jj);
                     }
         totmean_yXM(imp,ii) = totmean_yX[ii];
                }
     for (int ii=0;ii<NZ;ii++){
     for (int jj=0;jj<NZ;jj++){
         Sigma_W_yZM(ii+NZ*imp,jj) = Sigma_W_yZ(ii,jj);
         Sigma_B_yZM(ii+NZ*imp,jj) = Sigma_B_yZ(ii,jj);
                     }
         totmean_yZM(imp,ii) = totmean_yZ[ii];
                }

     arma::mat theta0 = res32["theta"];
     arma::mat Tmat0 = res32["Tmat"];
     arma::mat sig20 = res32["sig2"];
     iter2 = as<int>(res32["iter"]);
     iterM(gg,imp) = iter2;

     // Rcpp::Rcout << "a450" <<    std::flush << std::endl;

     //-------------- DATASET REPLICATE WEIGHTS --------------------------

     Rcpp::List res35=bifie_mla2_estimation_replicates( N__, NC__,
         wgttot__, wgtrep__, wgtlev1__, wgtlev2__,
         idcluster_table2, theta0, Tmat0, sig20, NX, NZ,  X__,
         Z__, y__,  globconv,  maxiter_rep, NP,
        recov_constraint, is_rcov_constraint, NRC );

     Rcpp::NumericMatrix pars2 = res35["pars_temp"]; // extract estimated parameters
     for (int pp=0;pp<NP;pp++){
         for (int rr=0;rr<RR;rr++){
            parsMrep( pp + gg*NP, rr+RR*imp ) = pars2(pp,rr);
                        }
                 }

     iter3 = res35["iter_temp"];
     for (int rr=0;rr<RR;rr++){
         iterMrep(gg,rr+RR*imp) = iter3[rr];
                 }

     // compute standard errors
     for (int pp=0;pp<NP;pp++){
        pars0[pp] = parsM( pp + gg*NP, imp );

        for (int rr=0;rr<RR;rr++){
               pars20(pp,rr) = parsMrep( pp + gg*NP, rr+RR*imp );
                           }
                    }

     Rcpp::List res41 = varjack_bias_helper( pars0, pars20, fayfac );
     Rcpp::NumericVector pars_var = res41["pars_var"];

     for (int pp=0;pp<NP;pp++){
         parsVar( pp + gg*NP, imp ) = pars_var[ pp ];
             }
     } // end loop groups
     //----------
     }
     /////////////////////// end imputations

     //////////////////////////////////////////////////////
     //*** post-processing

     ///*** Rubin inference
     Rcpp::List parsL = rubin_rules_univ( parsM, parsVar );
     Rcpp::Rcout << " " << std::endl;
     //*************************************************
     // OUTPUT
     return Rcpp::List::create(
             Rcpp::_["parsM"] = parsM,
             Rcpp::_["parsrepM"] = parsMrep,
             Rcpp::_["parsVar"] = parsVar,
             Rcpp::_["parsL"] = parsL,
                Rcpp::_["GG"] = GG,
             Rcpp::_["iterM"] = iterM,
              Rcpp::_["iterMrep"] = iterMrep,
             Rcpp::_["fvcovM"] = fvcovM,
              Rcpp::_["Npers"] = Npers,
             Rcpp::_["Nclusters"] = Nclusters,
              Rcpp::_["NP"] = NP,
             Rcpp::_["idcluster_table"] = idcluster_table2,
             Rcpp::_["Sigma_W_yXM"] = Sigma_W_yXM,
             Rcpp::_["Sigma_B_yXM"] = Sigma_B_yXM,
             Rcpp::_["Sigma_W_yZM"] = Sigma_W_yZM,
             Rcpp::_["Sigma_B_yZM"] = Sigma_B_yZM,
             Rcpp::_["totmean_yXM"] = totmean_yXM,
             Rcpp::_["totmean_yZM"] = totmean_yZM
              );
}

//*************************************************************************

//*************************************************************************
//  bifiesurvey_rcpp_replication_variance
// [[Rcpp::export]]
Rcpp::NumericVector bifiesurvey_rcpp_replication_variance( Rcpp::NumericVector pars,
        Rcpp::NumericMatrix pars_repl, Rcpp::NumericVector fay_factor )
{
    Rcpp::NumericVector res__ = varjack_helper( pars, pars_repl, fay_factor );
    return res__;
}
//*************************************************************************

//*************************************************************************
//  bifiesurvey_rcpp_rubin_rules
// [[Rcpp::export]]
Rcpp::List bifiesurvey_rcpp_rubin_rules( Rcpp::NumericMatrix estimates, Rcpp::NumericMatrix variances )
{
    Rcpp::List res__ = rubin_rules_univ( estimates, variances );
    return res__;
}
//*************************************************************************

// Rcout << "The value of chi2M : \n" << chi2M << "\n";
