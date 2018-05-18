//// File Name: univar_helpers.h
//// File Version: 7.22


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>
// #include <Rcpp.h>


using namespace Rcpp;
using namespace arma;



//***************************************
// path model
Rcpp::List bifie_pathmodel_helper( Rcpp::NumericMatrix dat1 ,
    Rcpp::NumericMatrix wgt , Rcpp::NumericVector group_values ,
    Rcpp::NumericVector group_index1 , Rcpp::NumericVector vars_index ,
    Rcpp::NumericMatrix R , int NP0 , Rcpp::NumericMatrix coeff_index ,
    Rcpp::NumericMatrix E , Rcpp::NumericMatrix L , int NL ,
    Rcpp::NumericVector L_row_index , Rcpp::NumericVector R_row_index ,
    Rcpp::NumericVector unreliability ){

    int N = dat1.nrow();
    int GG = group_values.size();
    int RR = wgt.ncol();
    int VV = vars_index.size();
    
    int group_index = group_index1[0] ;    
    Rcpp::NumericMatrix indcases(N,GG);
    Rcpp::NumericVector ncases(GG);
    Rcpp::NumericMatrix sumwgt1(GG,RR);
    int WW = wgt.ncol();
    int NR = R.nrow();

    // parameter vector
    int NP = 2*NP0 + 2*NR ;
    Rcpp::NumericMatrix parsM( GG*NP , WW ) ;
    int Npow = coeff_index.ncol();
    int Nrci = coeff_index.nrow();
    
    //**************************
    // extract usable cases
    for (int nn=0;nn<N;nn++){ // beg nn
       for (int gg=0; gg < GG ; gg++ ){  // beg gg
        if ( dat1(nn,group_index) == group_values[gg] ){ 
            indcases(nn,gg)=1 ;
                    }
        for (int vv=0;vv<VV;vv++){ // beg vv  independent vars
           if (  R_IsNA( dat1(nn, vars_index[vv] ) ) ){
            indcases(nn,gg)=0;
                    }
                } // end vv  independent vars
        if ( indcases(nn,gg) == 1 ){
            ncases[gg] ++ ;
            for (int ww=0;ww<WW; ww++){ 
                sumwgt1(gg,ww) += wgt(nn,ww) ;
                }
            break ;
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
        E_gg(_,vv) = E(_,vv) ;
                }
            
    for (int nn=0;nn<N;nn++){
    if ( indcases(nn,gg) == 1){
        for (int vv=0;vv<VV;vv++){
          if ( ! R_IsNA(vars_index[vv] ) ){        
              X(igg,vv) = dat1( nn , vars_index[vv] ) ;
                        }
                    }
        for (int ww=0;ww<WW;ww++){
            wgg(igg,ww) = wgt(nn,ww);
                }
        igg ++ ;
            }
         }

    // create derived variables
    if (NL>0){
        for (int jj=0;jj<NL;jj++){
        for (int nn=0;nn<ngg;nn++){
        for (int vv=0;vv<VV;vv++){
        if ( L(jj,vv) != 0 ){
            X( nn , L_row_index[jj] ) += X( nn , vv ) * L( jj , vv ) ; 
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
        weights_gg[gg] += wgg(nn,ww) ;
        for (int vv=0;vv<VV;vv++){
           means_gg(vv,gg) += wgg(nn,ww) * X(nn,vv) ;
           sds_gg(vv,gg) += wgg(nn,ww) * X(nn,vv) * X(nn,vv) ;
           for (int xx=vv+1;xx<VV;xx++){
              covs_gg( vv , xx+gg*VV) += wgg(nn,ww) * X(nn,vv) * X(nn,xx) ;  
                        }
                    }
                }        
    for (int vv=0;vv<VV;vv++){
        means_gg(vv,gg) = means_gg( vv,gg) / weights_gg[gg] ;
        sds_gg(vv,gg) = sqrt( ( sds_gg(vv,gg) - weights_gg[gg] * pow( means_gg( vv , gg) , 2.0 ) )/
                    weights_gg[gg] ) ; 
        covs_gg(vv,vv+gg*VV) = pow( sds_gg(vv,gg) , 2.0 ) ; 
                }
    
    for (int vv=0;vv<VV;vv++){
           for (int xx=vv+1;xx<VV;xx++){
             covs_gg( vv , xx+gg*VV) = covs_gg(vv,xx+gg*VV) - weights_gg[gg] * means_gg(vv,gg) * means_gg(xx,gg);
             covs_gg( vv , xx+gg*VV) = covs_gg(vv,xx+gg*VV) / weights_gg[gg] ;
             covs_gg( xx , vv+gg*VV) = covs_gg(vv,xx+gg*VV) ;
                        }
                }         
                    
    // calculate Cronbach's alpha for latent variables and adjust
    // measurement error variance
    for (int jj=0;jj<NL;jj++){
        NItems=0;
        for (int vv=0;vv<VV;vv++){
            vbar += sqrt( L(jj,vv) * L(jj,vv) ) * covs_gg( vv , vv+gg*VV) ;
            if ( L(jj,vv) != 0 ){
                NItems = NItems + 1 ;
                    }
            for (int xx=0;xx<VV;xx++){
               if ( xx != vv ){    
                cbar += L(jj,vv) * L(jj,xx) * covs_gg( vv , xx+gg*VV) ;
                        }            
                    }
                }
        vbar = vbar / NItems ;
        cbar = cbar / ( NItems*(NItems-1) ) ;
        alpha = NItems * cbar / ( vbar + ( NItems - 1 ) * cbar ) ;
        lli = L_row_index[jj];
        E_gg( lli , lli ) = ( 1 - alpha ) * covs_gg( lli , lli +gg*VV) ; 
            }
    
    // correction for unreliability of variances
    for (int vv=0;vv<VV;vv++){
        if ( unreliability[vv] > 0 ){
             E_gg( vv , vv ) = covs_gg( vv , vv + gg*VV ) * unreliability[vv] ;    
                    }
            }

            
    // adjustment of covariance matrix for measurement error covariance
    for (int vv=0;vv<VV;vv++){
      for (int xx=0;xx<VV;xx++){
          XtX_gg(vv,xx+gg*VV) = covs_gg( vv , xx+gg*VV) - E_gg(vv,xx);
                }
            }

    // calculate regression coefficients
    int nv_rr = 0 ;
    
    for (int rr=0;rr<NR;rr++){
        // int rr=2;
        
        // select predictors
        nv_rr=0 ;
        for (int vv=0;vv<VV;vv++){
            if ( ( R(rr,vv) != 0 ) ){
                preds_rr[nv_rr] = vv ;
                nv_rr = nv_rr+1 ; 
                    }        
                }
        
        arma::mat sigma_XX=arma::zeros(nv_rr,nv_rr);
        arma::mat sigma_XY=arma::zeros(nv_rr,1);
        arma::mat regr_coef_stand=arma::zeros(nv_rr,1);
                    
        for (int vv=0;vv < nv_rr; vv++){
        for (int xx=0;xx < nv_rr; xx++){    
            sigma_XX( vv, xx ) = XtX_gg( preds_rr[vv] , preds_rr[xx]+gg*VV) ;
                        }
            sigma_XY(vv,0) = XtX_gg( preds_rr[vv] , R_row_index[rr]+gg*VV) ;
                    }
                    
        arma::mat regr_coef = arma::mat( arma::pinv( sigma_XX ) * sigma_XY ) ;
        for (int vv=0 ; vv < nv_rr ; vv ++ ){
            regr_coef_stand(vv,0) = regr_coef(vv,0) * 
                sqrt( sigma_XX(vv,vv) ) /
                  sqrt( XtX_gg( R_row_index[rr], R_row_index[rr]+gg*VV) ) ;
            B_est_gg( R_row_index[rr] , preds_rr[vv] ) = regr_coef(vv,0);
            B_est_stand_gg( R_row_index[rr] , preds_rr[vv] ) = regr_coef_stand(vv,0);
                       }
        // calculate R-squared
        arma::mat var_expl = arma::mat( arma::trans( regr_coef ) * sigma_XX * regr_coef );               
        rsquared[rr] = var_expl(0,0) / XtX_gg( R_row_index[rr], R_row_index[rr]+gg*VV) ;
        residvar[rr] = XtX_gg( R_row_index[rr], R_row_index[rr]+gg*VV) * ( 1- rsquared[rr] );            
        }
            
    //--- calculate parameter list
    
    // extract coefficients
    for (int zz=0;zz<Nrci; zz++){ // beg zz
        parsM( zz + gg*NP , ww ) = 1;
        parsM( zz + NP0 + gg*NP , ww ) = 1;
        for (int hh=2;hh<Npow-1;hh++){
        if ( ! ( R_IsNA( coeff_index(zz,hh+1) ) ) ){
             parsM( zz + gg*NP , ww ) = parsM( zz + gg*NP , ww ) * 
                B_est_gg( coeff_index(zz,hh+1)-1  , coeff_index(zz,hh)-1 ) ;
             parsM( zz + NP0+ gg*NP , ww ) = parsM( zz + NP0 + gg*NP , ww ) * 
                B_est_stand_gg( coeff_index(zz,hh+1)-1  , coeff_index(zz,hh)-1 ) ;         
                            }
                    }
        if ( ! ( R_IsNA( coeff_index(zz , 0 ) ) ) ){
            parsM( coeff_index(zz , 0 ) - 1 + gg*NP  , ww ) += parsM( zz + gg*NP , ww ) ;
                parsM( coeff_index(zz , 0 ) - 1 + NP0 + gg*NP  , ww ) += parsM( zz + NP0+ gg*NP , ww ) ;            
                    }
        if ( ! ( R_IsNA( coeff_index(zz , 1 ) ) ) ){
            parsM( coeff_index(zz , 1 ) - 1 + gg*NP  , ww ) += parsM( zz + gg*NP , ww ) ;
                parsM( coeff_index(zz , 1 ) - 1 + NP0 + gg*NP  , ww ) += parsM( zz + NP0+ gg*NP , ww ) ;            
                    }                                                            
            }  // end zz
    // r squared
    for (int rr=0;rr<NR;rr++){
        parsM( rr + 2*NP0 + gg*NP , ww ) = rsquared[rr] ;
        parsM( rr + NR + 2*NP0 + gg*NP , ww ) = residvar[rr] ;
                }
            } // end ww
        }   // end gg
                                 
    //--- output     
    return Rcpp::List::create( 
        _["parsM"] = parsM ,
        _["ncases"] = ncases ,
        _["sumwgt1"] = sumwgt1
            ) ;
        }
//**********************************************************
// Rcpp::Rcout << "sigma_XX(0,0)= " <<  sigma_XX(0,0) << " " <<
//    "(1,1) " <<  sigma_XX(1,1)     <<  std::flush << std::endl ;        

//**********************************************************
// signum function
Rcpp::NumericVector bifie_sign( double x ){
    Rcpp::NumericVector y1(1);    
    if ( x > 0 ){ y1[0] = 1; }
    if ( x < 0 ){ y1[0] = -1 ; }
    return( y1 );
}
//**********************************************************

//**********************************************************
// converts a 1-column matrix into a vector
Rcpp::NumericVector matr2vec( Rcpp::NumericMatrix matr1){
    int N1=matr1.nrow();
    Rcpp::NumericVector vect1(N1) ;
    for (int zz=0;zz<N1;zz++){
        vect1[zz] = matr1(zz,0) ;
            }
    return(vect1)    ;
}
//**********************************************************

//**********************************************************
// matrix entries form parsM into a larger matrix
// starting from row zz in pars_full
Rcpp::List matrix_entry( Rcpp::NumericMatrix parsM , Rcpp::NumericMatrix pars_full1 ,
        int vv )
{
    int PP=parsM.nrow();
    int WW=parsM.ncol();    
    int zz=vv;    
    int N1 = pars_full1.nrow() ;
    Rcpp::NumericMatrix pars_full( N1, WW );
    pars_full = pars_full1 ;    
    for (int pp=0;pp<PP;pp++){
        for (int ww=0;ww<WW;ww++){    
           pars_full( zz,ww) = parsM(pp,ww) ;
                }
         zz++ ;
       }
    Rcpp::NumericVector zz2(1);
    zz2[0] = zz;
    return Rcpp::List::create( 
        _["pars_full"] = pars_full ,
        _["zz2"] = zz2 
            ) ;
}
//**********************************************************


//**********************************************************
// mean and SD in case of multiple groups
Rcpp::List univar_helper_multiple_V2group( Rcpp::NumericMatrix dat1, 
     Rcpp::NumericMatrix wgt1,  Rcpp::NumericVector vars_index ,
     Rcpp::NumericVector group_index1 , Rcpp::NumericVector group_values ){    

     int group_index = group_index1[0] ;
     int N = dat1.nrow() ;  
     int VV = vars_index.size() ;  
     int GG=group_values.size() ;
     int WW=wgt1.ncol() ;
     
     Rcpp::NumericMatrix sumwgt1(GG*VV,WW);  
     Rcpp::NumericMatrix ncases1(GG*VV,1);  
     double eps = 1E-20;
     
     // compute sum of weights  
     for (int nn=0;nn<N;nn++){ // beg nn
        for ( int gg=0;gg<GG;gg++){ // beg gg
            for (int vv=0;vv<VV;vv++){ // beg vv
                if ( ! R_IsNA( dat1(nn,vars_index[vv]) ) ){  // beg if IsNA                
                    if ( dat1(nn,group_index) == group_values[gg] ){ // beg if dat(,gg)==
                        for (int hh=0;hh<WW;hh++){  // beg hh        
                            sumwgt1(gg+vv*GG,hh) += wgt1(nn,hh) ;
                           } // end hh
                           ncases1(gg+vv*GG,0) ++ ;  
                           //               break ;                
                       } // end if dat(,gg) == group_values[gg]
                   } // end if IsNA
            } // end vv                  
           } // end gg
        }  // end nn
            
     int VV1=VV*GG ;       
     // compute means and standard deviations             
     Rcpp::NumericMatrix mean1(VV1,WW);  
     Rcpp::NumericMatrix sd1(VV1,WW);

     Rcpp::NumericMatrix mean1vv(GG,WW);  
     Rcpp::NumericMatrix sd1vv(GG,WW);  
       
     for (int vv=0; vv<VV;vv++){
         for (int hh=0;hh<WW;hh++){ 
             for (int gg=0;gg<GG;gg++){
                 mean1vv(gg,hh) = 0  ;  
                 sd1vv(gg,hh) = 0  ;
             }
         }
         for (int nn=0;nn<N;nn++){ // begin nn
             if ( ! R_IsNA( dat1(nn,vars_index[vv]) ) ){  // beg if IsNA                
                 for ( int gg=0;gg<GG;gg++){ // begin gg
                     if ( dat1(nn,group_index) == group_values[gg] ){     
                         for (int hh=0;hh<WW;hh++){        
                             mean1vv(gg,hh) += wgt1(nn,hh) * dat1(nn,vars_index[vv]) ;  
                             sd1vv(gg,hh) += wgt1(nn,hh) * pow( dat1(nn,vars_index[vv]) , 2.0) ;
                        }
                        break ;
                    } // end if    
                } // end gg
            } // end if IsNA                             
        }  // end nn
         for (int gg=0;gg<GG;gg++){  // begin gg
             for (int hh=0;hh<WW;hh++){ // beg hh            
                 mean1(vv*GG+gg,hh) = mean1vv(gg,hh) / ( sumwgt1(gg+vv*GG,hh) + eps );  
                 sd1(vv*GG+gg,hh) = sqrt( ( sd1vv(gg,hh) - sumwgt1(gg+vv*GG,hh)*
                     pow(mean1(vv*GG+gg,hh),2) ) /( sumwgt1(gg+vv*GG,hh) - 1 ) ) ;
            }
           } // end gg
           
     } // end vv
                        
    return Rcpp::List::create( 
        Rcpp::_["sumwgt1"] = sumwgt1 ,
        Rcpp::_["ncases1"] = ncases1 ,
        Rcpp::_["mean1"] = mean1 ,
        Rcpp::_["sd1"] = sd1
        ) ;                          
}    
//********************************************************************

     // Rcpp::Rcout << "c200 " << std::flush << std::endl ;

//**********************************************************
// compute standard errors using replicated statistics
Rcpp::NumericVector varjack_helper( Rcpp::NumericVector pars , 
    Rcpp::NumericMatrix pars_jack , Rcpp::NumericVector fayfac ){
    int PP=pars.size();
    // use int here in the subfunction
    int RR=pars_jack.ncol();
    Rcpp::NumericVector pars_var(PP) ;
    
    //@@fayfac 
    int NF = fayfac.size();
    double f1=0;
        //--
    
    double tmp1=0;    
    for (int pp=0; pp <PP ; pp++){
       tmp1=0;
       //@@fayfac
       f1 = fayfac[0];
       //--
       for (int rr=0;rr<RR;rr++){
          //@@fayfac
          if (NF>1){
             f1 = fayfac[rr];
                 }
          //--              
          tmp1 += f1 * pow( pars_jack(pp,rr) - pars[pp] ,2.0) ;
                }
       pars_var[pp] = tmp1 ;
            }
    return( pars_var ) ;
}
//**********************************************************

//**********************************************************
// compute standard errors and bias correction using replicated statistics
Rcpp::List varjack_bias_helper( Rcpp::NumericVector pars , 
    Rcpp::NumericMatrix pars_jack , Rcpp::NumericVector fayfac ){
    int PP=pars.size();
    // use int here in the subfunction
    int RR=pars_jack.ncol();
    Rcpp::NumericVector pars_var(PP) ;
    Rcpp::NumericVector pars_bias(PP) ;    

    //@@fayfac 
    int NF = fayfac.size();
    double f1=0;
        //--    
    
    double tmp1=0;    
    for (int pp=0; pp <PP ; pp++){
       tmp1=0;
       for (int rr=0;rr<RR;rr++){
          pars_bias[pp] += pars_jack(pp,rr) ;
                }       
       pars_bias[pp] = pars_bias[pp] / RR ;
       //@@fayfac
       f1 = fayfac[0];
       //--          
       for (int rr=0;rr<RR;rr++){
          //@@fayfac
          if (NF>1){
             f1 = fayfac[rr];
                 }              
          tmp1 += f1 * pow( pars_jack(pp,rr) - pars_bias[pp] ,2.0) ;
                }
       pars_var[pp] = tmp1 ;
            }
    return Rcpp::List::create( 
        _["pars_bias"] = pars_bias ,
        _["pars_var"] = pars_var      
                ) ;    
}        
        
//**********************************************************
// Rubin's rules for combining estimates
Rcpp::List rubin_rules_univ( Rcpp::NumericMatrix parsM , Rcpp::NumericMatrix pars_varM){

    int NP = parsM.nrow();
    int Nimp = parsM.ncol();    
    Rcpp::NumericVector pars(NP) ;
    Rcpp::NumericVector pars_varWithin(NP) ;
    Rcpp::NumericVector pars_varBetween(NP) ;
    Rcpp::NumericVector pars_se(NP) ;
    Rcpp::NumericVector pars_fmi(NP) ;
    
    double tmp1=0;
    double tmp2=0;
    double tmp3=0;
    double eps=1e-10;
    double Nimp2 = Nimp ;
    
    for ( int pp = 0; pp < NP ; pp++){
    tmp1=0;
    tmp2=0;  // squared variation , variance between
    tmp3=0; // variance within
    for ( int ii = 0; ii < Nimp ; ii++ ){
       tmp1 += parsM(pp,ii) ;
       tmp2 += parsM(pp,ii) * parsM(pp,ii) ;
       tmp3 += pars_varM(pp,ii) ;
            }
    
    pars[pp] = tmp1 / Nimp2 ;
    pars_varWithin[pp] = tmp3 / Nimp2 ;           
    pars_varBetween[pp] = ( tmp2 - Nimp2 * pow( pars[pp] , 2.0) ) / (Nimp2 - 1+eps ) ;
    // ARb 2014-09-10:  added "1.0" instead of "1"
    pars_se[pp] = sqrt( pars_varWithin[pp] + ( 1.0 + 1/Nimp2) * pars_varBetween[pp] ) ;
    pars_fmi[pp] = ( 1.0 + 1/Nimp2) * pars_varBetween[pp] / pow(pars_se[pp] + eps,2.0) ;
      }

return Rcpp::List::create( 
    _["pars"] = pars ,
    _["pars_se"] = pars_se ,     
    _["pars_varWithin"] = pars_varWithin ,
    _["pars_varBetween"] = pars_varBetween ,
    _["pars_fmi"] = pars_fmi 
            ) ;              
}
//**********************************************************

//**********************************************************
// subroutine frequency calculation
Rcpp::List bifiehelper_freq( Rcpp::NumericMatrix dat1 , Rcpp::NumericMatrix wgt ,
    Rcpp::NumericVector group_index1 , Rcpp::NumericVector group_values ,
    Rcpp::NumericVector vars_values_numb , Rcpp::NumericMatrix vars_values ,
    Rcpp::NumericVector vars_index , Rcpp::NumericVector vars_values_numb_cumsum ){

    int WW = wgt.ncol() ; 
    int N=dat1.nrow();
    int group_index = group_index1[0] ;
    int GG=group_values.size() ;
    int VV = vars_index.size() ;    
    int VV2= ( vars_values_numb_cumsum[VV] ) * GG ;

    Rcpp::NumericMatrix perc1(VV2,WW) ;
    Rcpp::NumericMatrix perc2(VV2,WW) ;
    Rcpp::NumericMatrix sumwgt(VV*GG,WW) ;
    Rcpp::NumericVector ncases(VV*GG) ;
    Rcpp::NumericVector ncases1(VV2) ;    
    
    int ind=0;    
    
    //******************************************************
    //********** count frequencies *************************
    for (int nn=0; nn<N ; nn++){ // beg nn
    for (int gg=0; gg<GG;gg++){  // beg gg
     if ( dat1(nn,group_index) == group_values[gg] ){  // beg if group_values
      for (int vv=0;vv<VV;vv++){   // beg vv       
      if (! R_IsNA(dat1(nn,vars_index[vv]) ) ){ // beg R_IsNA      
       for (int aa=0;aa < vars_values_numb[vv] ; aa++){  // beg aa
       if ( dat1(nn,vars_index[vv]) == vars_values(aa,vv) ){  // beg if vars_values
        ncases[vv*GG+gg] ++ ;          
        ind  = vars_values_numb_cumsum[vv]*GG + aa + gg * vars_values_numb[vv]  ;
        ncases1[ind] ++ ;   
        for (int ww=0;ww<WW;ww++){ // beg ww    
            perc1(ind,ww) += wgt(nn,ww) ;
            sumwgt(vv*GG+gg,ww) += wgt(nn,ww) ;
               }  // end ww
        break ;
            }  // end if vars_values                
          }   // end aa
      }  // end R_IsNA   
         } // end vv
       break ;
        }    // end if = group_values
      } // end gg            
    } // end nn
    
    //******************************************************
    //********** compute percentages ***********************
    for (int gg=0; gg<GG;gg++){  // beg gg
      for (int vv=0;vv<VV;vv++){   // beg vv         
       for (int aa=0;aa < vars_values_numb[vv] ; aa++){  // beg aa       
        ind  = vars_values_numb_cumsum[vv]*GG + aa + gg * vars_values_numb[vv]  ;
        for (int ww=0;ww<WW;ww++){ // beg ww    
            perc2(ind,ww) = perc1(ind,ww) / sumwgt(vv*GG+gg,ww) ; 
               }  // end ww            
          }   // end aa
         } // end vv
      } // end gg        
    
    
    return Rcpp::List::create( 
        _["sumwgt"] = sumwgt ,
        _["ncases1"] = ncases1 ,
        _["ncases"] = ncases  ,
        _["perc1"] = perc1 ,
        _["perc2"] = perc2   ,
        _["vars_values_numb_cumsum"] = vars_values_numb_cumsum
        ) ;  
}
//**********************************************************
// Rcpp::Rcout << "c200 " << std::flush << std::endl ;  


//**********************************************************
// correlation
Rcpp::List bifiehelpers_correl( Rcpp::NumericMatrix dat1 , Rcpp::NumericVector ind_cases ,
    Rcpp::NumericVector group_values , Rcpp::NumericVector group_index1 ,
    Rcpp::NumericMatrix wgt , Rcpp::NumericVector vars_index ,
    Rcpp::NumericMatrix itempair_index    
    ){

    int WW = wgt.ncol() ; 
    int N = wgt.nrow() ;
    int VV = vars_index.size() ;
    int group_index = group_index1[0] ;
    int GG=group_values.size() ;
    int ZZ = itempair_index.nrow() ;
    
    Rcpp::NumericMatrix mean1(VV*GG,WW) ;
    Rcpp::NumericMatrix sd1(VV*GG,WW) ;
    Rcpp::NumericMatrix sumwgt1(GG,WW) ;
    Rcpp::NumericVector ncases1(GG) ;
    
    // create index matrix for covariances and correlations
    Rcpp::NumericMatrix cov1(ZZ*GG,WW) ;
    Rcpp::NumericMatrix cor1(ZZ*GG,WW) ;
    
    
    for (int nn=0;nn<N ; nn++){  // beg nn
    if ( ind_cases[nn] == 1 ){  // beg ind_cases == 1     
    for (int gg=0;gg<GG;gg++){    
       if ( dat1(nn, group_index ) == group_values[gg] ){ // beg group == group_values[gg]
          ncases1[gg] ++ ;       
          for (int ww=0;ww<WW;ww++){  // beg ww       
        // weights  
        sumwgt1( gg ,ww ) += wgt(nn,ww) ;
        //**** loop variables means and SDs
        for (int vv=0;vv<VV;vv++){
           int tmpi2 = vv*GG + gg ;        
           mean1( tmpi2 , ww ) += wgt(nn,ww) * dat1( nn , vars_index[vv] ) ;
           sd1( tmpi2 , ww ) += wgt(nn,ww) * pow(dat1( nn , vars_index[vv] ),2.0) ;    
                    }
        
        // covariances
        for (int zz=0;zz<ZZ;zz++){ // beg zz item pairs
           int vv1=itempair_index(zz,0) ;
           int vv2=itempair_index(zz,1) ;
           if (vv1<vv2){
             int tmpi1 = zz*GG + gg ;               
             cov1(tmpi1 , ww ) += wgt(nn,ww) * dat1( nn , vars_index[vv1] ) * 
                  dat1( nn , vars_index[vv2] ) ;
                 }
                } // end zz (item pairs)
                 } // end ww
        break ;     
                } // end group == group_values[gg]            
            } // end gg
              } // end if ind_cases == 1
            }   // end nn
            
    // calculate descriptive statistics
    for (int ww=0;ww<WW;ww++){
    for (int vv=0;vv<VV; vv++){
       for (int gg=0;gg<GG;gg++){
           int zz=vv*GG+gg ;       
           mean1( zz , ww ) = mean1( zz , ww ) / sumwgt1( gg , ww ) ;           
           sd1(zz, ww) = sqrt( ( sd1(zz,ww) - sumwgt1( gg , ww ) * pow(mean1(zz,ww) ,2.0) )/
                 ( sumwgt1( gg , ww ) - 1 ) ) ;
                    }
                    }       
                    
    for (int pp=0;pp<ZZ;pp++){
        if (  itempair_index(pp,0) == itempair_index(pp,1) ){
          int ipp = itempair_index(pp,0) ;        
          for (int gg=0;gg<GG;gg++){
        cov1( pp*GG+gg , ww ) = pow( sd1( ipp*GG+gg , ww) , 2.0)  ;                   
                    }
                }
            }
    
            
    // calculate covariance and correlation
    for (int zz=0;zz<ZZ;zz++){
       int vv1=itempair_index(zz,0) ;
       int vv2=itempair_index(zz,1) ;
       if (vv1 < vv2 ){        
       for (int gg=0;gg<GG;gg++){
           int tmpi1 = zz*GG + gg ;         
        cov1(tmpi1,ww) = ( cov1(tmpi1,ww) - sumwgt1(gg,ww) * mean1( vv1*GG + gg , ww) * 
            mean1( vv2*GG + gg , ww ) ) / ( sumwgt1( gg , ww ) - 1 )  ;
                }
              }
            }
    
    // calculate correlations                
    for (int zz=0 ; zz < ZZ;zz++){
      int vv1=itempair_index(zz,0) ;
      int vv2=itempair_index(zz,1) ;
      for (int gg=0; gg<GG ;gg++){
          cor1(zz*GG+gg,ww) = cov1(zz*GG+gg,ww) / sd1( vv1*GG+gg,ww) / sd1( vv2*GG+gg,ww) ;
            }
        }
    }
        
    //**** output        
    return Rcpp::List::create( 
        _["sumwgt1"] = sumwgt1 , 
        _["mean1"] = mean1 ,
        _["sd1"] = sd1 , 
        _["cov1"] = cov1 ,
        _["cor1"] = cor1 ,
        _["ncases1"] = ncases1
        ) ;  
}
//**********************************************************

//**********************************************************
// linear regression
Rcpp::List bifiehelpers_linreg( Rcpp::NumericMatrix dat1 , 
    Rcpp::NumericVector group_values , Rcpp::NumericVector dep_index ,
    Rcpp::NumericVector pre_index , Rcpp::NumericMatrix wgt ,
    Rcpp::NumericVector group_index1 ){
        
int WW = wgt.ncol() ;
int N = wgt.nrow() ;
int VV = pre_index.size() ;
int group_index = group_index1[0] ;
int GG=group_values.size() ;    
int VV2=(2*VV+2)*GG;

Rcpp::NumericMatrix sumwgt1(GG,WW);
Rcpp::NumericVector ncases(GG);
Rcpp::NumericMatrix regr_coef(VV2,WW) ;
Rcpp::NumericMatrix indcases(N,GG);

//**************************
// extract usable cases
for (int nn=0;nn<N;nn++){ // beg nn
   for (int gg=0; gg < GG ; gg++ ){  // beg gg
    if ( dat1(nn,group_index) == group_values[gg] ){ 
        indcases(nn,gg)=1 ;
                }
    if (  R_IsNA( dat1(nn, dep_index[0] ) ) ){ // beg R_IsNA dep var
        indcases(nn,gg)=0;
                } // end R_IsNA dep var
    for (int vv=0;vv<VV;vv++){ // beg vv  independent vars
       if (  R_IsNA( dat1(nn, pre_index[vv] ) ) ){
        indcases(nn,gg)=0;
                }
            } // end vv  independent vars
        if ( indcases(nn,gg) == 1 ){
            ncases[gg] ++ ;
            for (int ww=0;ww<WW; ww++){ 
                sumwgt1(gg,ww) += wgt(nn,ww) ;
                    }
            break ;
                }
        }   // end gg
  }   // end nn    


double sig2=0;
double sig2a=0;

for (int gg=0;gg<GG;gg++){
    int igg=0;
     
    // create design matrices
    int ngg = ncases[gg] ;
    arma::mat X0=arma::zeros( ngg , VV ) ;
    arma::mat X=arma::zeros( ngg , VV ) ;
    arma::colvec y0=arma::zeros( ngg ) ;
    arma::colvec y=arma::zeros( ngg ) ;
    arma::mat wgt0=arma::zeros( ngg , WW) ;
    
    Rcpp::NumericMatrix M_pre(VV,WW) ;
    Rcpp::NumericMatrix SD_pre(VV,WW) ;
    Rcpp::NumericMatrix M_dep(1,WW) ;
    Rcpp::NumericMatrix SD_dep(1,WW) ;
    
    //***** define input matrices
    for (int nn=0;nn < N ; nn++){  // beg nn
        if ( indcases(nn,gg)==1){  // beg indcases gg 
            y0(igg,0) = dat1(nn, dep_index[0] ) ;
            for (int vv=0;vv<VV;vv++){  // beg vv
               X0(igg,vv) = dat1(nn,pre_index[vv] ) ;
                    }  // end vv                        
            igg ++ ;
            } // end if indcases gg
          } // end nn
    
    //**** define used matrices
    double wtmp=0;
    
    for ( int ww=0; ww <WW; ww++){ // beg ww
    igg=0;
    for (int nn=0;nn < N ; nn++){  // beg nn
        if ( indcases(nn,gg)==1){  // beg indcases gg
            wgt0(igg,ww) = wgt(nn,ww);        
            wtmp = sqrt( wgt(nn,ww)) ;    
            y(igg,0) = y0(igg,0)*wtmp ;
            M_dep(0,ww) += y0(igg,0) * wgt(nn,ww) ;
            SD_dep(0,ww) += y0(igg,0) * y0(igg,0) * wgt(nn,ww)  ;        
            for (int vv=0;vv<VV;vv++){  // beg vv
               X(igg,vv) = X0(igg, vv )*wtmp ;
               M_pre(vv,ww) += X0(igg,vv) * wgt(nn,ww) ;
               SD_pre(vv,ww) += X0(igg,vv) * X0(igg,vv) * wgt(nn,ww) ;               
                    }  // end vv
            igg ++ ;        
            } // end if indcases gg
          } // end nn      
    
    //*** fit linear model      
    arma::colvec coef = arma::solve(X, y);      // fit model y ~ X
    // arma::colvec resid = y - X*coef;            // residuals
    // not that the weights are already included in the residual calculation
    arma::colvec resid(ngg);
    for (int hh=0;hh<ngg;hh++){
        resid(hh,0) = y0(hh,0) ;
        for (int vv=0;vv<VV;vv++){
           resid(hh,0) = resid(hh,0) - X0(hh,vv) * coef(vv,0) ;
            }    
    }
         
    // sig2 = arma::as_scalar( arma::trans(resid)*resid );
    sig2=0;
    for (int hh=0;hh<ngg;hh++){
          sig2 += pow( resid(hh,0),2.0) * wgt0(hh,ww) ;    
       }
    double sggww = sumwgt1(gg,ww) ;
    sig2a = sig2 / ( sggww - VV ) ;
    
    // collect all regression coefficients
    // unstandardized coefficients
    for (int vv=0;vv<VV;vv++){
        regr_coef( vv + gg*(2*VV+2) , ww ) = coef(vv,0) ;
            }
    regr_coef( VV + gg*(2*VV+2) , ww ) = sqrt( sig2a ) ;  // sigma          
    // compute R^2
    M_dep(0,ww) = M_dep(0,ww) / sggww ;
    double sig3 = SD_dep(0,ww) - sggww * pow( M_dep(0,ww) , 2.0) ; 
    SD_dep(0,ww) = sqrt( sig3 / ( sggww - 1 ) ) ;
    regr_coef( VV+1 + gg*(2*VV+2) , ww ) = 1 - sig2 / sig3 ;
    
    // compute standardized coefficients
    for (int vv=0;vv<VV;vv++){
        M_pre(vv,ww) = M_pre(vv,ww) / sggww ;
        SD_pre(vv,ww) = SD_pre(vv,ww) - sggww *  pow( M_pre(vv,ww) , 2.0) ;
        SD_pre(vv,ww) = sqrt( SD_pre(vv,ww) / ( sggww - 1 ) ) ;
        regr_coef( VV+2+vv + gg*(2*VV+2) , ww ) = 
           coef(vv,0) / SD_dep(0,ww) * SD_pre(vv,ww) ;
                }
    } // end ww
} // end gg
       
return Rcpp::List::create( 
    _["ncases"] = ncases ,
    _["sumwgt1"] = sumwgt1 ,
    _["regr_coef"] = regr_coef 
    ) ; 
}
//**********************************************************


//**********************************************************
// BIFIE helpers Wald test
Rcpp::List bifiehelpers_waldtest( int VV , Rcpp::NumericVector Ccols ,
    Rcpp::NumericMatrix parsM , Rcpp::NumericMatrix parsrepM ,
    int ii , int RR , Rcpp::NumericVector fayfac ,
    arma::mat ACdes , arma::colvec Ardes ){
    
    //@@fayfac 
    int NF = fayfac.size();
    double f1=0;
        //--

    //*** calculate covariance matrix for imputation ii
    arma::mat var_w = arma::zeros(VV,VV);
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
       var_w(vv1,vv2) += f1 * ( parsrepM( Ccols[vv1] , rr+ii*RR ) - parsM( Ccols[vv1] , ii ) )*
            ( parsrepM( Ccols[vv2] , rr+ii*RR ) - parsM( Ccols[vv2] , ii ) ) ;
                }    
    // var_w(vv1,vv2) = var_w(vv1,vv2) ;
    var_w(vv2,vv1) = var_w(vv1,vv2) ;
                }
            }
    
    //***  compute chi squared Wald test statistic
    // compute covariance matrix of hypothesis
    arma::mat var_hyp = arma::mat( ACdes * var_w * trans( ACdes) ) ;
        
    // compute inverse of variance matrix of hypothesis
    arma::mat var_hypinv = arma::inv( var_hyp ) ;
    // parameter vector
    arma::colvec parm_vec= arma::zeros(VV,1);
    for (int vv=0;vv<VV;vv++){
        parm_vec(vv,0) = parsM( Ccols[vv] , ii ) ;
                }
    // hypothesis statistic
    arma::mat hyp_stat = arma::mat( ACdes * parm_vec - Ardes ) ;
    arma::mat chi2 = arma::mat( trans( hyp_stat ) * var_hypinv * hyp_stat ) ;             
    
    return Rcpp::List::create( 
        _["chi2"] = chi2 ,
        _["var_w"] = var_w
        ) ; 
}
//**********************************************************
    

//**********************************************************
// BIFIE helpers Wald test
Rcpp::List bifiehelpers_waldtest_vcov( int VV , Rcpp::NumericVector Ccols ,
    Rcpp::NumericMatrix parsM , Rcpp::NumericMatrix parsrepM ,
    int ii , int RR , Rcpp::NumericVector fayfac ,
    arma::mat ACdes , arma::colvec Ardes ){
    
    //@@fayfac 
    int NF = fayfac.size();
    double f1=0;
        //--c

    //*** calculate covariance matrix for imputation ii
    arma::mat var_w = arma::zeros(VV,VV);
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
       var_w(vv1,vv2) += f1 * ( parsrepM( Ccols[vv1] , rr+ii*RR ) - parsM( Ccols[vv1] , ii ) )*
            ( parsrepM( Ccols[vv2] , rr+ii*RR ) - parsM( Ccols[vv2] , ii ) ) ;
                }    
//    var_w(vv1,vv2) = fayfac[0] * var_w(vv1,vv2) ;
    var_w(vv2,vv1) = var_w(vv1,vv2) ;
                }
            }
    
    //***  compute chi squared Wald test statistic
    // compute covariance matrix of hypothesis
    arma::mat var_hyp = arma::mat( ACdes * var_w * trans( ACdes) ) ;
    
    // compute inverse of variance matrix of hypothesis
    arma::mat var_hypinv = arma::mat( var_hyp ) ;
    // parameter vector
    arma::colvec parm_vec= arma::zeros(VV,1);
    for (int vv=0;vv<VV;vv++){
        parm_vec(vv,0) = parsM( Ccols[vv] , ii ) ;
                }
    // hypothesis statistic
    arma::mat hyp_stat = arma::mat( ACdes * parm_vec - Ardes ) ;
    arma::mat chi2 = arma::mat( trans( hyp_stat ) * var_hypinv * hyp_stat ) ;             
    
    // int df=ACdes.n_rows;        
    return Rcpp::List::create( 
        _["chi2"] = chi2 ,
        _["var_w"] = var_w ,
        _["var_hyp"] = var_hyp ,
        _["hyp_stat"] = hyp_stat
        ) ; 
}
//**********************************************************    
    
//**********************************************************
// etasquared and d statistics
Rcpp::List bifiehelpers_etasquared( Rcpp::NumericMatrix mean1M ,
    Rcpp::NumericMatrix sd1M , Rcpp::NumericMatrix sumweightM , int GG ){
    
    int WW = sd1M.ncol() ;
    // calculate total mean
    Rcpp::NumericVector totmean(WW) ; 
    Rcpp::NumericVector sumwgt(WW) ; 
    Rcpp::NumericVector expl_var(WW) ;
    Rcpp::NumericVector resid_var(WW) ;
    Rcpp::NumericMatrix eta2(1,WW) ;
    int GG2 = GG * (GG - 1 ) / 2 ;
    Rcpp::NumericMatrix dstat(GG2 , WW ) ;
    

    for ( int ww=0; ww < WW ; ww++){ // beg ww
        for (int gg=0;gg<GG; gg++){ // beg gg
          sumwgt[ww] += sumweightM(gg,ww) ;
          totmean[ww] += sumweightM(gg,ww) * mean1M(gg,ww) ;      
// Rcpp::Rcout << "mean1M(gg,0)=" << mean1M(gg,ww) << std::flush << std::endl ;

            }  // end gg
        totmean[ww] = totmean[ww] / sumwgt[ww] ;  
        for (int gg=0;gg<GG; gg++){  // beg gg
          expl_var[ww] += sumweightM(gg,ww)*pow( mean1M(gg,ww) - totmean[ww] , 2.0 ) ;
          resid_var[ww] += (sumweightM(gg,ww)-1)*pow( sd1M(gg,ww) , 2.0 ) ; 
          eta2(0,ww) = sqrt( expl_var[ww] / ( expl_var[ww] + resid_var[ww] ) ) ;
            }  // end gg
        
        // calculate d statistics
        int ii=0;
        for ( int gg1=0; gg1 < GG - 1 ; gg1++){
        for (int gg2=gg1+1; gg2 < GG ; gg2++){
           dstat(ii,ww) = mean1M(gg1,ww) - mean1M(gg2,ww) ;
           dstat(ii,ww) = dstat(ii,ww) / sqrt( 0.5 * ( pow(sd1M(gg1,ww),2.0) + 
                    pow(sd1M(gg2,ww),2.0) ) ) ;
           ii++ ;
                }
            }
     }  // end ww


    return Rcpp::List::create( 
        _["eta2"] = eta2 ,
        _["dstat"] = dstat
        ) ;   
}
//**********************************************************  
  
//**********************************************************
// cross tabulation
Rcpp::List bifiehelpers_crosstab( Rcpp::NumericMatrix dat1 , Rcpp::NumericMatrix wgt ,
    Rcpp::NumericVector group_values , Rcpp::NumericVector group_index1 , 
    Rcpp::NumericVector vars_values1 , Rcpp::NumericVector vars_index1 ,
    Rcpp::NumericVector vars_values2 , Rcpp::NumericVector vars_index2 ,    
    Rcpp::NumericMatrix design_pars 
        ){
   int N=dat1.nrow();   
   int group_index = group_index1[0] ;        
   int GG=group_values.size() ;
   int WW=wgt.ncol();
   int VV1 = vars_values1.size() ;
   int VV2 = vars_values2.size() ;
   int ZZ = VV1*VV2*GG ;
   int zz=0;
   
   Rcpp::NumericMatrix ncases( ZZ , 1 ) ;
   Rcpp::NumericMatrix sumwgt( ZZ , WW ) ;
    for (int nn=0;nn<N; nn++){  // beg nn
    for (int gg=0;gg<GG;gg++){ // beg gg
    if ( dat1( nn , group_index ) == group_values[gg] ){  // beg if groupval gg
    for (int vv1=0;vv1<VV1;vv1++){   // beg vv1
       for (int vv2=0;vv2<VV2;vv2++){    // beg vv2
          if ( ( dat1(nn, vars_index1[0]) == vars_values1[vv1] ) &
               ( dat1(nn, vars_index2[0]) == vars_values2[vv2] ) ){  // beg if comb
            ncases( vv2 + vv1*VV2 + gg*VV1*VV2 , 0 ) ++ ;
            for (int ww=0;ww<WW;ww++){  // beg ww
               sumwgt( vv2 + vv1*VV2 + gg*VV1*VV2 , ww ) += wgt(nn,ww) ; 
                  } // end ww
                        break ;                   
                }  // end if combination of values
            }  // end vv2
        }   // end vv1
        } // end if groupval gg
    } // end gg
    }  // end nn
    
    Rcpp::NumericMatrix ncases_gg(GG,1);
    Rcpp::NumericMatrix sumwgt_gg(GG,WW);
    
    Rcpp::NumericMatrix probs_joint(ZZ,WW) ;
    Rcpp::NumericMatrix probs_rowcond(ZZ,WW) ; 
    Rcpp::NumericMatrix probs_colcond(ZZ,WW) ;
    
    Rcpp::NumericMatrix probs_rowmarg(GG*VV1,WW);
    Rcpp::NumericMatrix probs_colmarg(GG*VV2,WW);
    
    //**** compute counts within each group value gg
    for (int gg=0;gg<GG;gg++){ // beg gg
    for (int vv1=0;vv1<VV1;vv1++){  // beg vv1
       for (int vv2=0;vv2<VV2;vv2++){  // beg vv2
           ncases_gg(gg,0) += ncases( vv2 + vv1*VV2 + gg*VV1*VV2 , 0 ) ;
          for (int ww=0;ww<WW;ww++){
            sumwgt_gg(gg,ww) += sumwgt( vv2 + vv1*VV2 + gg*VV1*VV2 , ww ) ;
                    }
                    }  // end vv2
                }  // end vv1
           } // end gg
    
    //**** compute joint probabilities
    int ind=0;
    for (int gg=0;gg<GG;gg++){ // beg gg
    for (int vv1=0;vv1<VV1;vv1++){  // beg vv1
       for (int vv2=0;vv2<VV2;vv2++){  // beg vv2
          ind = vv2 + vv1*VV2 + gg*VV1*VV2 ; 
          for (int ww=0;ww<WW;ww++){
             probs_joint(ind,ww)= sumwgt( ind , ww ) / sumwgt_gg( gg ,ww);
                    }
                }  // end vv2
           }  // end vv1
        } // end gg       
    //**** compute conditional rowwise probabilities          
    double tmp1=0;
    for (int ww=0;ww<WW;ww++){
    for ( int gg=0;gg<GG;gg++){
    for (int vv1=0;vv1<VV1;vv1++){ // beg vv1
        tmp1=0;
        for (int vv2=0;vv2<VV2;vv2++){
           ind = vv2 + vv1*VV2 + gg*VV1*VV2 ;
           tmp1 += sumwgt( ind , ww ) ;
                    }
           probs_rowmarg( vv1 + gg*VV1 , ww ) = tmp1 / sumwgt_gg( gg , ww ) ;  
        for (int vv2=0;vv2<VV2;vv2++){
           ind = vv2 + vv1*VV2 + gg*VV1*VV2 ;
           probs_rowcond( ind , ww ) = sumwgt( ind , ww ) / tmp1 ;
                    }
        } // end vv1    
    } //end gg
    } // end ww
    
    //**** compute conditional columnwise probabilities
    for (int ww=0;ww<WW;ww++){
    for (int gg=0;gg<GG;gg++){
    for (int vv2=0;vv2<VV2;vv2++){
    tmp1=0;
         for (int vv1=0;vv1<VV1;vv1++){
        ind = vv2 + vv1*VV2 + gg*VV1*VV2 ;
        tmp1 += sumwgt( ind , ww ) ;
                } // end vv1
         probs_colmarg( vv2 + gg*VV2 , ww ) = tmp1 / sumwgt_gg(gg,ww) ;
         for (int vv1=0;vv1<VV1;vv1++){
           ind = vv2 + vv1*VV2 + gg*VV1*VV2 ;
           probs_colcond( ind , ww ) = sumwgt( ind , ww ) / tmp1 ;
                    }            
    }   // end vv2
    } // end gg
    } // end ww
           
    //*** effect size w  (also include Cramer's V)
    Rcpp::NumericMatrix w_es(2*GG,WW);  // w and Cramer's V
    int VV3 = VV1 ;
    if (VV2 < VV1 ){ VV3 = VV2 ; }
    double tmp2=0;
    for (int ww=0; ww<WW;ww++){ // beg ww
    for (int gg=0;gg<GG;gg++){ // beg gg
    for (int vv1=0;vv1<VV1;vv1++){ // beg vv1
    for (int vv2=0;vv2<VV2;vv2++){  // beg vv2
       ind = vv2 + vv1*VV2 + gg*VV1*VV2 ;
       // probability under independence;
       tmp2 = probs_rowmarg( vv1 + gg*VV1 , ww ) * probs_colmarg( vv2 + gg*VV2, ww );
       w_es(gg,ww) += pow( probs_joint( ind , ww ) - tmp2 , 2.0 ) / tmp2 ;
    //  Rcpp::Rcout << "w_es " <<  w_es(gg,0) << "  gg,ww" << gg << "  " << ww <<  std::flush << std::endl ;   
                    } // end vv2
         } // end vv1
    w_es(gg,ww) = sqrt( w_es(gg,ww) ) ;  // effect size w
    w_es(GG+gg,ww ) = w_es(gg,ww) / sqrt( VV3 - 1.0 ) ; // Cramer's V
    } // end gg
    } // end ww
    
    //*** gamma measure of association (Goodman)
    Rcpp::NumericMatrix gamma_es( GG , WW ) ;
    
    double pc=0;
    double pd=0;
    // double h1=0;
    double h2=0;
    int VV12 = VV1*VV2 ;
    
    for (int ww=0;ww<WW;ww++){   // beg ww
    for (int gg=0; gg < GG ; gg++){ // beg gg
        pc=0;
        pd=0;
        // numer := 0
        //for i:=2..N do
        //    for j:=1..(i-1) do
        //        numer := numer + sign(x[i] - x[j]) * sign(y[i] - y[j])
        // return numer
        
        for (int zzi = 0 ; zzi < VV12 ; zzi++ ){ // beg zzi
        for( int zzj = 0 ; zzj < VV12 ; zzj++ ){ // beg zzj
        
        double h1 = vars_values1[ design_pars(zzi , 3) ] - vars_values1[ design_pars(zzj , 3) ] ;
        Rcpp::NumericVector h1s = bifie_sign( h1 ) ;
        
        h2 = vars_values2[ design_pars(zzi , 4) ] > vars_values2[ design_pars(zzj , 4) ] ; 
        Rcpp::NumericVector h2s = bifie_sign( h2 ) ;
        
        double h3 = h1s[0] * h2s[0] ;
        if ( zzi== zzj ){ h3 = 0 ; }
        int vv1 = design_pars(zzj , 3) ;
        int vv2 = design_pars(zzj , 4) ;
        ind = vv2 + vv1*VV2 + gg*VV1*VV2 ;   
        
        double p3= probs_joint( design_pars(zzi,4) + design_pars(zzi , 3)*VV2 + gg*VV1*VV2 , ww) ;        
                
        // concordant cases
        if ( h3 > .99 ){
            pc += probs_joint(ind,ww) * p3;
           }  // end if concordant
        
        // discordant cases
        if ( h3 < -.99 ){ 
            pd += probs_joint(ind,ww) * p3 ;
           }  // end if discordant   
        }  // end zzj
        } // end zzi
        gamma_es(gg,ww) = ( pc - pd ) / ( pc + pd ) ;
    } // end gg
    } // end ww
    //**** lambda PRE measures
    Rcpp::NumericMatrix lambda(3*GG, WW) ;  // lambda_V1, lambda_V2, lambda
    double l1y=0;
    double l2y=0;
    double p1=0;
    double t11=0;
    double l1x=0;
    double l2x=0;
    
    for (int ww=0; ww<WW;ww++){ // beg ww
    for (int gg=0;gg<GG;gg++){ // beg gg
    
        //-- lambda_Y
        l1y=0;
        l2y=0;
        for (int vv1 = 0 ;vv1<VV1;vv1++){
            t11 = 0 ;
            for (int vv2=0;vv2<VV2;vv2++){
                 ind = vv2 + vv1*VV2 + gg*VV1*VV2 ; 
                 p1 = probs_joint( ind , ww ) ;    
                 if (p1> t11 ){ t11 = p1 ; } 
                }
            l1y += t11 ;
            }    
        for (int vv2=0;vv2<VV2;vv2++){
           p1 = probs_colmarg( vv2 + gg*VV2 , ww ) ;
           if (p1 > l2y ){ l2y = p1 ; }
                    }
        lambda( 2 + gg*3 , ww ) = ( l1y - l2y ) / ( 1 - l2y ) ; 
        
        //-- lambda_X
        l1x=0;
        l2x=0;
        for (int vv2 = 0 ;vv2<VV2;vv2++){
            t11 = 0 ;
            for (int vv1=0;vv1<VV1;vv1++){
                 ind = vv2 + vv1*VV2 + gg*VV1*VV2 ; 
                 p1 = probs_joint( ind , ww ) ;    
                 if (p1> t11 ){ t11 = p1 ; } 
                }
            l1x += t11 ;
            }    
        for (int vv1=0;vv1<VV1;vv1++){
           p1 = probs_rowmarg( vv1 + gg*VV1 , ww ) ;
           if (p1 > l2x ){ l2x = p1 ; }                
           } 
        lambda( 1 + gg*3 , ww ) = ( l1x - l2x ) / ( 1 - l2x ) ; 
        
        //-- symmetric lambda measure
        lambda( 0 + gg*3 , ww ) = ( l1x + l1y - l2x - l2y ) / ( 2 - l2x - l2y ) ;
    
        
    } // end gg
    } // end ww
    
    //**** Kruskal's tau    
    Rcpp::NumericMatrix kruskal_tau(3*GG,WW) ;
    
    double t1y=0 ;
    double t2y=0 ;
    double t1x=0;
    double t2x=0;
    
    
    for (int ww=0;ww<WW;ww++){ // beg ww
    for (int gg=0;gg<GG;gg++){ // beg gg
    //*** tau_Y
    t1y=0;
    t2y=0;
    for (int vv2=0;vv2<VV2;vv2++){  // beg vv1
        t2y += pow( probs_colmarg( vv2 + gg*VV2 , ww ) , 2.0 );
                  }  // end vv2
    for ( int vv1=0; vv1 < VV1 ; vv1++){  // beg vv1
    for (int vv2=0;vv2<VV2;vv2++){   // beg vv2
       ind = vv2 + vv1*VV2 + gg*VV1*VV2 ; 
       t1y += pow( probs_joint( ind , ww ) , 2.0 ) / probs_rowmarg( vv1+gg*VV1, ww ) ;
                }   // end vv2
            }  // end vv1
    kruskal_tau(2+gg*3,ww) = ( t1y - t2y ) / ( 1 - t2y ) ;
    
    //*** tau_X
    t1x=0;
    t2x=0;
    for (int vv1=0;vv1<VV1;vv1++){
        t2x += pow( probs_rowmarg( vv1 + gg*VV1 , ww ) , 2.0 );
                  }
    for ( int vv2=0; vv2 < VV2 ; vv2++){  // beg vv2
    for (int vv1=0;vv1<VV1;vv1++){   // beg vv1
       ind = vv2 + vv1*VV2 + gg*VV1*VV2 ; 
       t1x += pow( probs_joint( ind , ww ) , 2.0 ) / probs_colmarg( vv2+gg*VV2, ww ) ;
                }   // end vv1
            }  // end vv2
    kruskal_tau(1+gg*3,ww) = ( t1x - t2x ) / ( 1 - t2x ) ;
    
    //*** symmetric tau
    kruskal_tau(0+gg*3,ww) = ( t1x+t1y-t2x-t2y)/(2-t2x-t2y);
    }  // end gg
    } // end ww
    
    //**** output vector of parameters from crosstable
    
    // probs_joint    ZZ
    // probs_rowcond  ZZ
    // probs_colcond   ZZ
    // probs_rowmarg   VV1*GG
    // probs_colmarg  VV2*GG
    // w_es  2*GG
    // gamma_es  GG
    // lambda   3*GG
    // kruskal_tau  3*GG
    
    int CTP = 3*ZZ + VV1*GG + VV2*GG + 2*GG + GG + 3*GG + 3*GG ;
    Rcpp::NumericMatrix crosstab_pars( CTP , WW) ;
    zz=0;
    //*** probs_joint
    Rcpp::List res2= matrix_entry( probs_joint , crosstab_pars ,zz ) ;   
    Rcpp::NumericVector zz2 = res2["zz2"] ;
    zz=zz2[0] ;
    Rcpp::NumericMatrix v2 = res2["pars_full"] ;
    //*** rowcond
    res2= matrix_entry( probs_rowcond , v2 ,zz ) ;
    zz2 = res2["zz2"] ;
    zz=zz2[0] ;
    Rcpp::NumericMatrix v3 = res2["pars_full"] ;
    //*** colcond
    res2= matrix_entry( probs_colcond , v3 ,zz ) ;
    zz2 = res2["zz2"] ;
    zz=zz2[0] ;
    Rcpp::NumericMatrix v4 = res2["pars_full"] ;
    //*** rowmarg
    res2= matrix_entry( probs_rowmarg , v4 ,zz ) ;
    zz2 = res2["zz2"] ;
    zz=zz2[0] ;
    Rcpp::NumericMatrix v5 = res2["pars_full"] ;
    //*** colmarg
    res2= matrix_entry( probs_colmarg , v5 ,zz ) ;
    zz2 = res2["zz2"] ;
    zz=zz2[0] ;
    Rcpp::NumericMatrix v6 = res2["pars_full"] ;
    //*** w_es
    res2= matrix_entry( w_es , v6 ,zz ) ;
    zz2 = res2["zz2"] ;
    zz=zz2[0] ;
    Rcpp::NumericMatrix v7 = res2["pars_full"] ;
    //*** gamma_es
    res2= matrix_entry( gamma_es , v7 ,zz ) ;
    zz2 = res2["zz2"] ;
    zz=zz2[0] ;
    Rcpp::NumericMatrix v8 = res2["pars_full"] ;
    //*** lambda
    res2= matrix_entry( lambda , v8 ,zz ) ;
    zz2 = res2["zz2"] ;
    zz=zz2[0] ;
    Rcpp::NumericMatrix v9 = res2["pars_full"] ;
    //*** kruskal_tau
    res2= matrix_entry( kruskal_tau , v9 ,zz ) ;
    zz2 = res2["zz2"] ;
    zz=zz2[0] ;
    Rcpp::NumericMatrix v10 = res2["pars_full"] ;    
    //*************************************************    
    // OUTPUT                    
    return Rcpp::List::create( 
        _["ncases"] = ncases ,
        _["ncases_gg"] = ncases_gg , 
        _["sumwgt"] = sumwgt ,
        _["sumwgt_gg"] = sumwgt_gg ,
        _["crosstab_pars"] = v10 
        ) ;  
}
//**********************************************************

    
//**********************************************************    
// bifie_helper_ecdf
Rcpp::NumericVector bifie_helper_ecdf( Rcpp::NumericMatrix dat1 , 
    Rcpp::NumericVector wgt , Rcpp::NumericVector breaks ,
    Rcpp::NumericVector group_values , Rcpp::NumericVector group_index1 ,
    Rcpp::NumericVector vars_index , int ii ,
    Rcpp::NumericMatrix ncasesM , Rcpp::NumericMatrix sumwgtM ,
    int maxval , int quanttype ){
        
    int N=dat1.nrow();    
    int GG=group_values.size() ;
    int BB=breaks.size() ;
    int VV=vars_index.size();
    Rcpp::NumericVector ecdf_temp(BB);
    arma::colvec vals_temp(N) ;
    arma::colvec wgt_temp(N) ;
    arma::colvec group_temp(N) ;
    arma::uvec indvv(N);
    arma::colvec wgt_tempsort(N);
    arma::colvec vals_tempsort(N);
    arma::colvec group_tempsort(N);
    arma::colvec wgt_tempsort_gg(N);
    arma::colvec vals_tempsort_gg(N);
    int group_index = group_index1[0] ;
    int ZZ=VV*BB*GG;
    Rcpp::NumericVector ecdfMtemp(ZZ);
    
    Rcpp::NumericMatrix a1(N,4) ;
    double eps = 1e-20 ;
    double maxval2 = maxval ;
    
    int uu=0;
    int ngg=0;
    int hh=0;
    int MM=0;
    
    //*** count observed responses
    for (int nn=0;nn<N;nn++){ // beg nn
    for (int gg=0;gg<GG;gg++){ // beg gg
    if ( dat1(nn,group_index ) == group_values[gg] ){ // beg val gg
    for (int vv=0;vv<VV;vv++){ // beg vv
        if ( ! R_IsNA( dat1(nn,vars_index[vv] ) ) ){    
        ncasesM(gg + vv*GG , ii ) ++ ;
        sumwgtM(gg+vv*GG,ii) += wgt(nn,0) ;
                    }
                }  // end vv
        }  // enf if val gg
    } // end gg
    } // end nn


    int bb=0;
    
    for (int vv=0;vv<VV;vv++){
    
    for (int nn=0;nn<N;nn++){ // beg nn
       group_temp(nn,0) = -1 ;
       wgt_temp(nn,0) = 0 ;
       vals_temp(nn,0) = maxval2 ;
       for (int gg=0;gg<GG;gg++){
       if ( dat1(nn,group_index ) == group_values[gg] ){ // beg val gg
          if ( ! R_IsNA( dat1(nn,vars_index[vv] ) ) ){ // beg non NA
               vals_temp(nn,0) = dat1(nn,vars_index[vv] ) ;
               wgt_temp(nn,0) = wgt(nn,0 ) ; 
               group_temp(nn,0) = gg ;
                } // end non NA
         } // end val gg
       } // end gg
    } // end nn 
                    
    // sort values
    indvv = arma::sort_index( vals_temp ) ;
    vals_tempsort =vals_temp.rows( indvv ) ;
    wgt_tempsort =wgt_temp.rows( indvv ) ;
    group_tempsort =group_temp.rows( indvv ) ;
    
    for ( int gg=0;gg<GG;gg++){
    ngg = ncasesM( gg + vv*GG , ii ) ;
    uu=0;
    for (int nn=0;nn<N;nn++){ // beg nn
       if ( group_tempsort(nn,0) == gg ){    // beg group = gg
         vals_tempsort_gg(uu,0) = vals_tempsort(nn,0) ;
         wgt_tempsort_gg(uu,0) = wgt_tempsort(nn,0) / sumwgtM( gg+vv*GG, ii) ;
         uu ++ ;
          } // end if group = gg
    }  // end nn
    
    
    if (ngg<N){        // beg if ngg < N
        for (int uu=ngg; uu<N;uu++){  // beg uu
         vals_tempsort_gg(uu,0) = 0 ;
         wgt_tempsort_gg(uu,0) = 0 ;
            } // end uu
        } // end if ngg < N
        
    // define ECDF    
    a1(0,0)=0;
    a1(0,1)=wgt_tempsort_gg(0,0);
    a1(0,2)=vals_tempsort_gg(0,0);
    a1(0,3)=vals_tempsort_gg(0,0);
    
    for (int uuu=1;uuu<ngg;uuu++){
        a1(uuu,0) = a1(uuu-1,1) ;
        a1(uuu,1) = a1(uuu-1,1) + wgt_tempsort_gg(uuu,0) ;
        a1(uuu,2) = a1(uuu-1,3);
        a1(uuu,3) = vals_tempsort_gg(uuu,0) ;
    }                
    uu=0;        
    bb=0;
    MM=0;
    while ( ( bb < BB ) & ( ! ( uu >= ngg ) ) ){
        if ( a1(uu,1) > breaks[bb] ){
           if (quanttype==1){    
            ecdf_temp[bb] = a1(uu,2) + ( a1(uu,3) - a1(uu,2) ) * 
                        ( breaks[bb] - a1(uu,0)  ) / ( a1(uu,1) - a1(uu,0) + eps) ;
            MM=bb;
                    }        
           if (quanttype==2){    
            ecdf_temp[bb] = a1(uu,2) ;
            MM=bb;
                    }                
            bb ++ ;
                    }  else {
            uu ++ ;
                    }
            } // end while
    // Rcpp::Rcout << "MM= " <<  MM <<  std::flush << std::endl ;        
       if ( MM < BB  ){
          hh=MM+1 ;
          while (hh < BB ){               
           ecdf_temp[bb] = a1(ngg-1,3) ;
           hh ++ ;
            }
           }
    
        for (int bb=0;bb<BB;bb++){
                ecdfMtemp[ bb + gg*BB + vv*GG*BB ] = ecdf_temp[bb] ;        
            }  // end bb
            
    }  // end gg
    
    } // end vv
    
    return ( ecdfMtemp ) ;
}
//**********************************************************    
    
//**********************************************************
//**** logistic regression    
Rcpp::List bifie_estlogistic_helper( Rcpp::NumericVector y ,
    Rcpp::NumericMatrix X, Rcpp::NumericVector wgt ,
    Rcpp::NumericVector beta0 , double eps , int maxiter ){
  int N=X.nrow();
  int P=X.ncol();
  double t1=0;
  
  //*** create matrices in Armadillo
  // design matrix X
  arma::mat Xa0(N,P);
  arma::mat Xa(N,P);
  for (int nn=0;nn<N;nn++){
     for (int pp=0;pp<P;pp++){
      Xa0(nn,pp) = X(nn,pp) ;
              }
           }
  // outcome matrix y
  arma::colvec ya(N) ;
  for (int nn=0;nn<N;nn++){
      ya(nn,0) = y[nn] ;
           }
  // regression coefficients
  arma::colvec beta_old(P);
  arma::colvec beta_new(P);
  for (int pp=0;pp<P;pp++){
      beta_old(pp,0) = beta0[pp] ;
              }
  
  // temporary values in iterations
  arma::colvec pred_logit(N);
  arma::colvec ypred(N);
  arma::colvec z(N);
  arma::colvec AM(N);
  arma::colvec wgta(N);
  double pardiff=100;
  
  int ii=0;
  
  while( ( pardiff > eps ) & ( ii < maxiter ) ){
  
  // within an iteration
  
  // calculate predicted logit value and probability
  for( int nn=0; nn <N; nn++){
  pred_logit(nn,0)=0;
  for ( int pp=0; pp <P; pp++){
      pred_logit(nn,0) += Xa0(nn,pp) * beta_old(pp,0) ;
          }
  if ( pred_logit(nn,0) < - 15 ){
      pred_logit(nn,0) = - 15 ;
                  }    
  ypred(nn,0) = 1 / ( 1 + exp( - pred_logit(nn,0) ) ) ;
  }
  // calculate entries for A matrix and outcome z
  for (int nn=0;nn<N;nn++){
      AM(nn,0) = ypred(nn,0) * ( 1 - ypred(nn,0) ) ;
      wgta(nn,0) = sqrt( AM(nn,0) * wgt[nn] ) ;
      z(nn,0) = pred_logit(nn,0) + ( ya(nn,0) - ypred(nn,0) )/AM(nn,0) ;
      z(nn,0) = wgta(nn,0) * z(nn,0) ;
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
      t1 = beta_old(pp,0) - beta_new(pp,0) ;
      if (t1 < 0 ){ t1 = -t1 ; }
      if (t1 > pardiff){ pardiff = t1 ; }
          }
  for (int pp=0;pp<P; pp++){
      beta_old(pp,0) = beta_new(pp,0);
          }
     ii ++ ;        
      }
      
  //*****************    
  // OUTPUT                    
  return Rcpp::List::create( 
      _["pardiff"] = pardiff ,
      _["beta"] = beta_new ,
      _["iter"] = ii
      ) ;  
}
//**********************************************************




//**********************************************************
// convert Rcpp matrix into Armadillo matrix
// !!!!!! change output to matrix format
Rcpp::List rcppmat2armamat( Rcpp::NumericMatrix matr ){
    int N = matr.nrow();
    int NX = matr.ncol() ;
    arma::mat Xa(matr.begin(), N, NX, false);     
    return Rcpp::List::create( 
       _["armamat"] = Xa
        ) ;  
}
//**********************************************************

//**********************************************************
// standard errors multilevel models
Rcpp::List mla2_se_fixed( arma::mat theta , arma::mat Tmat , arma::mat sig2 ,
        arma::mat Var_theta_rj , arma::mat AfjArj , 
        Rcpp::NumericVector wgtlev2 , arma::mat Xt2 ){

    int NX = theta.n_rows ;
    int NZ = Tmat.n_rows;
    int NC = Var_theta_rj.n_rows / NZ ;     
    arma::mat fvcov = arma::zeros(NX,NX);
    arma::mat mat2 = arma::zeros(NX,NX);
    arma::mat mat3 = arma::zeros(NX,NZ);
    arma::mat Cjinv = arma::zeros(NZ,NZ);
    arma::mat mat4 = arma::zeros(NZ,NX);
    arma::mat mat30 = arma::zeros(NX,NX);
    
    for (int jj=0;jj<NC;jj++){  // beg cluster jj
        // Var_theta_rj( ii1 + jj*NZ , ii2) = sig2(0,0) * Cjinv(ii1,ii2);
        for (int ii=0;ii<NZ;ii++){
        for (int hh=0;hh<NZ;hh++){
           Cjinv(ii,hh) = wgtlev2[jj] * Var_theta_rj( ii + jj*NZ , hh ) / sig2(0,0) ;
                    }
                }
        for (int ii=0;ii<NX;ii++){
        for (int hh=0;hh<NZ;hh++){
            mat3(ii,hh) = AfjArj( ii + jj*NX , hh ) ;
                    }
                }
        mat4 = arma::trans(mat3) ;
        mat30 = arma::mat( mat3 * Cjinv * mat4 ) ;
        mat2 = arma::mat( mat2 + mat30 ) ;
            }  // end cluster jj
    
    fvcov = arma::inv( arma::mat( Xt2 - mat2 ) );
    for (int ii=0;ii<NX;ii++){
    for (int hh=0;hh<NX;hh++){
        fvcov(ii,hh) = sig2(0,0) * fvcov(ii,hh) ;
                }
            }
    return Rcpp::List::create( 
       _["fvcov"] = fvcov
        ) ;  
}    
//**********************************************************


//**********************************************************
// multilevel within-between decomposition
Rcpp::List mla2_decomp( Rcpp::NumericMatrix V , 
    Rcpp::NumericMatrix idcluster_table , Rcpp::NumericVector wgttot_ ){
    
    int N = V.nrow();
    int NV = V.ncol();
    int NC = idcluster_table.nrow();
    Rcpp::NumericMatrix grmean(NC,NV);
    Rcpp::NumericVector grwgt(NC);
    
    Rcpp::NumericVector totmean(NV);
    Rcpp::NumericMatrix MSA(NV,NV);
    Rcpp::NumericMatrix MSE(NV,NV);
    Rcpp::NumericMatrix Sigma_W(NV,NV);
    Rcpp::NumericMatrix Sigma_B(NV,NV);     
    
    Rcpp::NumericVector wgttot(N);
    double w0=0;
    
    for (int nn=0;nn<N;nn++){
        w0 += wgttot_[nn] ;
                }
    for (int nn=0;nn<N;nn++){
        wgttot[nn] = wgttot_[nn] * N / w0 ;
                }
        

    int Ntot = 0 ; 
    int Gtot = 0 ;
    int sumng2 = 0 ;
    double h1=0;
    double eps2 = 1E-10 ;
    
    for (int jj=0;jj<NC;jj++){
        for (int nn=idcluster_table(jj,0); nn < idcluster_table(jj,1)+1 ; nn++){
            grwgt[jj] += wgttot[nn] ;
            for (int ii=0;ii<NV;ii++){
            h1 = wgttot[nn] * V(nn,ii) ;    
            grmean(jj,ii) += h1 ;
            totmean[ii] += h1 ;
                }
            }
        for (int ii=0;ii<NV;ii++){
            grmean(jj,ii) = grmean(jj,ii) / ( grwgt[jj] + eps2 );
                }
        Ntot += grwgt[jj] ;
        sumng2 += pow( grwgt[jj] , 2.0 ) ;    
        
            }                
            
    Gtot = NC * Ntot / N ;
        
    //------
    // calculate total mean
    for (int ii=0; ii <NV ; ii++ ){
       totmean[ii] = totmean[ii] / Ntot ;
            }
    //----
    // calculate variance components SSA
    for (int ii=0;ii<NV; ii++){
    for (int hh=ii; hh < NV ; hh++ ){
        for (int jj=0;jj<NC;jj++){
           MSA(ii,hh) += grwgt[jj] * ( grmean(jj,ii) - totmean[ii] ) *( grmean(jj,hh) - totmean[hh] ) ;           
                    }
    //       MSA(ii,hh) = MSA(ii,hh) / ( Ntot - Gtot ) ;
           MSA(ii,hh) = MSA(ii,hh) / ( Gtot - 1) ;
        if ( hh > ii ){
            MSA( hh , ii ) = MSA( ii , hh ) ;
            }
            }
         }

    //---
    // calculate variance components SSE
    
    double m1=0;
    
    for (int ii=0;ii<NV;ii++){
    for (int hh=ii;hh<NV;hh++){
        
        for (int jj=0;jj<NC;jj++){
        for (int nn= idcluster_table(jj,0); nn < idcluster_table(jj,1) +1 ; nn++){
          m1 = wgttot[nn] * ( V(nn,ii) - grmean(jj,ii) ) * ( V(nn,hh) - grmean(jj,hh) ) ;
          MSE(ii,hh) += wgttot[nn] * m1 ;
                }
            }
        // MSE(ii,hh) = MSE(ii,hh) / ( Ntot - 1 ) ;
        MSE(ii,hh) = MSE(ii,hh) / ( Ntot - Gtot ) ;
        if (hh>ii){
            MSE(hh,ii) = MSE(ii,hh);
               }     
        }
         }
    //--- estimation variance components
        double eps0 = Ntot * 1E-10 ;    
    for (int ii=0;ii<NV;ii++){
    for (int hh=0;hh<NV;hh++){
        Sigma_W(ii,hh) = MSE(ii,hh) ;
        Sigma_B(ii,hh) = ( MSA(ii,hh) - MSE(ii,hh) ) *
         ( Gtot - 1.0 ) / ( Ntot - sumng2 / ( Ntot + eps0 ) ) ;
        if ( ( MSA(hh,hh) < eps0 ) | ( MSA(ii,ii) < eps0 ) ){ 
                Sigma_B(ii,hh) = 0 ;
                                }
                }
            }    
 // Rcpp::Rcout << "MSA(0,0)=" << MSA(0,0) <<   std::flush << std::endl ;            
    return Rcpp::List::create( 
       _["totmean"] = totmean , _["Sigma_B"] = Sigma_B ,
       _["Sigma_W"] = Sigma_W ,
       _["grmean"] = grmean , _["grwgt"] = grwgt 
        ) ;      
}
//**********************************************************


//**********************************************************
// maximal difference of two armadillo matrices
Rcpp::NumericVector maxabsval_arma( arma::mat Tmat , arma::mat Tmat0 ){
    
    int N1 = Tmat.n_rows;
    int N2 = Tmat.n_cols;
    
    double absval=0;
    double temp=0;
    
    for (int rr=0;rr<N1;rr++){
    for (int cc=0;cc<N2;cc++){
        temp= Tmat(rr,cc) - Tmat0(rr,cc) ;
        if ( temp < 0 ){ 
            temp = - temp ;
                }
        if ( temp > absval ){
            absval = temp ;
                }
            }
        }
    Rcpp::NumericVector absval2(1);
    absval2[0] = absval ;
    return ( absval2  );
}
//**********************************************************    
    
//**********************************************************  
// compute deviation between succeeding iterations
Rcpp::NumericVector mla2_checkconv( arma::mat theta , arma::mat theta0 ,
    arma::mat Tmat , arma::mat Tmat0 , arma::mat sig2 ,
    arma::mat sig20 ){
    //--
    arma::mat absval2 = maxabsval_arma( theta , theta0 ) ;
    double absval=0;
    absval = absval2[0] ;
    absval2 = maxabsval_arma( Tmat , Tmat0 ) ;
    if (absval2[0] > absval ){
        absval = absval2[0] ;
                }    
    absval2 = maxabsval_arma( sig2 , sig20 ) ;
    if (absval2[0] > absval ){
        absval = absval2[0] ;
                }
    Rcpp::NumericVector absval3(1);
    absval3[0] = absval ;
    return(absval3 ) ;                
}                
//**********************************************************  

//**********************************************************
// variance decomposition Snijders and Bosker (2012)
Rcpp::NumericVector mla2_vardec( arma::mat theta , arma::mat Tmat , arma::mat sig2 ,
    Rcpp::NumericMatrix Sigma_B_yX ,Rcpp::NumericMatrix Sigma_W_yX ,
    Rcpp::NumericMatrix Sigma_B_yZ ,Rcpp::NumericMatrix Sigma_W_yZ ,
    Rcpp::NumericVector totmean_yZ ){
    
        int NX = theta.n_rows;
        int NZ = Tmat.n_rows;    

    Rcpp::NumericVector vardec(13);
    int ee=0;
    
    //** fixed effects Level 2
    ee=0; // index
    for (int ii=1;ii<NX;ii++ ){
    for (int hh=1;hh<NX;hh++){
       vardec[ee] += theta(ii,0) * theta(hh,0) * Sigma_B_yX( ii , hh ) ; 
                }
            }
    //** fixed effects Level 1
    ee=3;
    for (int ii=1;ii<NX;ii++ ){
    for (int hh=1;hh<NX;hh++){
       vardec[ee] += theta(ii,0) * theta(hh,0) * Sigma_W_yX( ii , hh ) ; 
                }
            }
    
    //** random effects Level 2
    ee=1;
    for (int ii=1;ii<NZ;ii++){
    for (int kk=1;kk<NZ;kk++){
       vardec[ee] += Tmat(ii,kk) * Sigma_B_yZ( kk , ii ) ;
                }
            }
    //** random effects Level 1
    ee=4;
    for (int ii=1;ii<NZ;ii++){
    for (int kk=1;kk<NZ;kk++){
       vardec[ee] += Tmat(ii,kk) * Sigma_W_yZ( kk , ii ) ;
                }
            }
    
    //** unexplained variance Level 2
    ee=2;
    vardec[ee] += Tmat(0,0);
    for (int kk=1;kk<NZ;kk++){
       // vardec[ee] += totmean_yZ[kk] * Tmat(kk,0);
       vardec[ee] += 2 * totmean_yZ[kk] * Tmat(kk,0);
                }
    for (int kk=1;kk<NZ;kk++){
    for (int hh=1;hh<NZ;hh++){    
       vardec[ee] += totmean_yZ[kk] * totmean_yZ[hh] * Tmat(kk,hh);
                }
            }
    
    //** unexplained variance Level 1
    ee=5;
    vardec[ee] = sig2(0,0);
    //** total variance
    ee=6;
    for (int hh=0;hh<ee;hh++){
        vardec[ee] += vardec[hh] ;
                }
    //** R2 Level 2
    ee=7;
    vardec[ee] = ( vardec[0] + vardec[1] ) / ( vardec[0] + vardec[1] + vardec[2] ) ;
    //** R2 Level 1
    ee=8;
    vardec[ee] = ( vardec[3] + vardec[4] ) / ( vardec[3] + vardec[4] + vardec[5] ) ;
    //** R2 total
    ee=9;
    vardec[ee] = 1 - (vardec[2] + vardec[5] ) / vardec[6] ;
    //** ICC unconditional
    ee=10;
    vardec[ee] = ( vardec[0] + vardec[1] + vardec[2] ) / vardec[6] ;
    //** ICC unconditional within-between decomposition
    ee=11;    
        vardec[ee] = Sigma_B_yZ(0,0) / ( Sigma_B_yZ(0,0) + Sigma_W_yZ(0,0) );
    //** ICC conditional
    ee=12;
    vardec[ee] = ( vardec[2] ) / ( vardec[2] + vardec[5] ) ;
        
    return vardec;    
    }
//**********************************************************    


//**********************************************************
// post-processing variance decompositions
Rcpp::List mla2_postproc( int N , int NX , int NZ , Rcpp::NumericVector y ,
    Rcpp::NumericMatrix X , Rcpp::NumericMatrix Z , 
    Rcpp::NumericMatrix idcluster_table ,
    Rcpp::NumericVector wgttot , arma::mat theta , arma::mat Tmat ,
    arma::mat sig2 , arma::mat Var_theta_rj , arma::mat  AfjArj , 
    Rcpp::NumericVector wgtlev2 , arma::mat Xt2){    
    
    // covariance decomposition of y and X
    Rcpp::NumericMatrix V(N,NX);
    for (int nn=0;nn<N;nn++){
        V(nn,0) = y[nn] ;
      if (NX>1){
         for (int ii=1;ii<NX;ii++){
         V(nn,ii) = X(nn,ii) ; 
                }
            }
        }
    Rcpp::List res11 = mla2_decomp( V , idcluster_table , wgttot ) ;                     
    Rcpp::NumericVector totmean_yX = res11["totmean"] ;
    Rcpp::NumericMatrix Sigma_B_yX = res11["Sigma_B"] ;
    Rcpp::NumericMatrix Sigma_W_yX = res11["Sigma_W"] ;
    
    //--- covariance decomposition of y and Z
    Rcpp::NumericMatrix V1(N,NZ);
    for (int nn=0;nn<N;nn++){
        V1(nn,0) = y[nn] ;
     if (NZ>1){
      for (int ii=1;ii<NZ;ii++){
         V1(nn,ii) = Z(nn,ii) ;
                }
            }
                }
                
    res11 = mla2_decomp( V1 , idcluster_table , wgttot ) ;                     
    Rcpp::NumericVector totmean_yZ = res11["totmean"] ;
    Rcpp::NumericMatrix Sigma_B_yZ = res11["Sigma_B"] ;
    Rcpp::NumericMatrix Sigma_W_yZ = res11["Sigma_W"] ;

    //--- variance decomposition
    Rcpp::NumericVector vardec = mla2_vardec( theta , Tmat , sig2 ,
        Sigma_B_yX ,Sigma_W_yX , Sigma_B_yZ , Sigma_W_yZ , totmean_yZ ) ;    
    
    arma::mat fvcov = mla2_se_fixed( theta , Tmat , sig2 , Var_theta_rj , 
                   AfjArj , wgtlev2 , Xt2 )["fvcov"];     
    
    return Rcpp::List::create( 
       Rcpp::_["totmean_yX"] = totmean_yX , 
       Rcpp::_["Sigma_B_yX"] = Sigma_B_yX ,
       _["Sigma_W_yX"] = Sigma_W_yX ,
       _["totmean_yZ"] = totmean_yZ , 
       _["Sigma_B_yZ"] = Sigma_B_yZ ,
       _["Sigma_W_yZ"] = Sigma_W_yZ  ,
       _["vardec"] = vardec , 
       _["fvcov"] = fvcov
                    ) ;  
        
    }
//**********************************************************    
    
//**********************************************************
// dummy variables
Rcpp::List create_dummies_mla2( int GG , Rcpp::NumericVector group ,
    Rcpp::NumericMatrix X , Rcpp::NumericMatrix Z , Rcpp::NumericVector y )
{ 
    int NX = X.ncol();
    int NZ = Z.ncol();
    int N = group.size() ;
    Rcpp::NumericMatrix dummy_inds(N,GG);
    Rcpp::NumericVector N_group(GG);

    int grnn=0;
    for (int nn=0;nn<N;nn++){
        grnn = group[nn] ;
        dummy_inds( nn , grnn ) = 1 ;
        if ( R_IsNA( y[nn] ) ){
            dummy_inds(nn,grnn) = 0 ;
        }
        for (int zz=0;zz<NX;zz++){
            if ( R_IsNA( X(nn,zz) ) ){
                dummy_inds(nn,grnn) = 0 ;
            }
        }
        for (int zz=0;zz<NZ;zz++){
            if ( R_IsNA( Z(nn,zz) ) ){
                dummy_inds(nn,grnn) = 0 ;
            }
        }
    }
    for (int nn=0;nn<N;nn++){
        for (int gg=0;gg<GG;gg++){
            N_group[ gg ] += dummy_inds(nn,gg) ;
        }
    }
    //--- output
    return Rcpp::List::create( 
            Rcpp::Named("N_group") = N_group ,  
            Rcpp::Named("dummy_inds") = dummy_inds 
            );
}
//**********************************************************    

//**********************************************************
// create cluster table
Rcpp::NumericMatrix create_idclustertable( Rcpp::NumericVector group ,
    Rcpp::NumericVector cluster , int NC)
{
    int N = group.size();
    Rcpp::NumericMatrix idcluster_table2(NC,4);    
    double cl=-99999;
    int jj=-1;    
    for (int nn=0;nn<N;nn++){
        if ( cluster[nn] > cl ){
            jj ++ ;
            cl = cluster[nn] ;
            idcluster_table2(jj,0) = nn ;
            idcluster_table2(jj,2) = cluster[nn] ;
            idcluster_table2(jj,3) = group[nn] ;
        } else {
            idcluster_table2(jj,1) = nn ;
        }
    }
    return idcluster_table2;
}
//**********************************************************        


//**********************************************************
// rescale level 1 weights
Rcpp::NumericVector rescale_lev1weights( Rcpp::NumericMatrix idcluster_table ,
        Rcpp::NumericVector wgtlev1 )
{
    int NC = idcluster_table.nrow();
    int N = wgtlev1.size();
    Rcpp::NumericMatrix wgtlev1_table(NC,2);    
    Rcpp::NumericVector wgtlev1a(N);    
    for (int jj=0;jj<NC;jj++){
        for (int ii=idcluster_table(jj,0) ; ii <idcluster_table(jj,1)+1 ; ii++){
            wgtlev1_table(jj,0) ++ ;
            wgtlev1_table(jj,1) += wgtlev1[ii] ;
        }
        for (int ii=idcluster_table(jj,0) ; ii <idcluster_table(jj,1)+1 ; ii++){
            wgtlev1a[ii] = wgtlev1[ii] / wgtlev1_table(jj,1) * wgtlev1_table(jj,0) ;     
        }
    }
     return wgtlev1a;
}
//**********************************************************


//**********************************************************
// E and M step of the algorithm
Rcpp::List mla2_emsteps( Rcpp::NumericMatrix X, Rcpp::NumericMatrix Z ,
    Rcpp::NumericMatrix idcluster_table , arma::mat Tmat ,
    arma::mat Xa , arma::mat theta , Rcpp::NumericVector y,
    Rcpp::NumericVector wgtlev2 , Rcpp::NumericVector wgtlev1 ,
    arma::mat ArjArj , arma::mat AfjArj , arma::mat sig2 ,
    Rcpp::NumericVector W1, Rcpp::NumericVector W2 ,
    arma::mat XtX , arma::mat Xty  ){

    int NZ = Z.ncol();                                                   
    int NX = X.ncol();
    int N = X.nrow();
    int NC = idcluster_table.nrow();    
        
    //--- EM iterations
    arma::mat Eterm1 = arma::zeros(NX,1);
    arma::mat Eterm2 = arma::zeros(NZ,NZ);
    arma::mat Eterm3 = arma::zeros(1,1);
    arma::mat theta_rj = arma::zeros(NC,NZ);
    arma::mat Var_theta_rj = arma::zeros(NC*NZ,NZ);

    double eps = 1E-10 ;
    // eps = 0 ;
    
    arma::mat Tmat0(NZ,NZ);
    
    for (int zz=0;zz<NZ;zz++){
    for (int hh=0;hh<NZ;hh++){    
        Tmat0(zz,hh) = Tmat(zz,hh) ;
        if (zz==hh){
           Tmat0(zz,hh) = Tmat0(zz,hh) + eps ;
                  }    
            }
        }
    
    arma::mat Tmatinv = arma::inv(Tmat0) ;
    
    arma::mat Cj = arma::zeros(NZ,NZ);
    arma::mat Cjinv = arma::zeros(NZ,NZ);
    arma::mat trj = arma::zeros(NZ,1);
    arma::mat theta_rj_jj = arma::zeros(NZ,1);
    arma::mat A1 = arma::zeros(NZ,NZ);
    double sumtr3 = 0 ;
    double rj=0;
    
    arma::mat ypred_ranef = arma::zeros(N,1);    
    // prediction fixed effects
    arma::mat ypred = arma::mat( Xa * theta );
        
    for (int jj=0;jj<NC;jj++){
    
        // id.jj <- idcluster_table[jj,]
        // id.jj <- seq( id.jj[1] , id.jj[2] )
        //  # compute Cj
        //  Cj <- ArjArj[ 1:NZ + (jj-1)*NZ , ] + sig2 * Tmatinv
        //        Cjinv <- solve(Cj)
        
        for (int ii1=0;ii1<NZ;ii1++){
        for (int ii2=0;ii2<NZ;ii2++){
            Cj(ii1,ii2) = ArjArj( ii1 + jj * NZ , ii2 ) +  sig2(0,0) * Tmatinv(ii1,ii2) ;
            if (ii1 == ii2 ){
                Cj(ii1,ii1) = Cj(ii1,ii1) + eps ;
                      }
                    }
                }
        Cjinv = arma::pinv(Cj) ;
        
        // # compute theta_rj^\ast
        // r1 <- y[ id.jj ] - X[ id.jj , ] %*% theta
        // theta_rj.jj <- Cjinv %*% t( Z[ id.jj , ] ) %*% wgtlev1M[id.jj,id.jj] %*% r1    
        // theta_rj[jj,] <- as.vector( theta_rj.jj )
        // Var_theta_rj[ 1:NZ + (jj-1)*NZ , ] <- sig2 * Cjinv        
        
        for (int ii=0;ii<NZ;ii++){
        trj(ii,0) = 0 ;
        for (int nn = idcluster_table(jj,0);nn < idcluster_table(jj,1)+1;nn++){
            trj(ii,0) += Z( nn , ii ) * wgtlev1[nn] * ( y[nn] - ypred(nn,0) ) ;
                }
            }
        theta_rj_jj = arma::mat( Cjinv * trj );
        for (int zz=0;zz<NZ;zz++){
            theta_rj(jj,zz) = theta_rj_jj(zz,0) ;
                     }
                     
        // Var_theta_rj[ 1:NZ + (jj-1)*NZ , ] <- sig2 * Cjinv             
        
        for (int ii1=0;ii1<NZ;ii1++){
        for (int ii2=0;ii2<NZ;ii2++){
           Var_theta_rj( ii1 + jj*NZ , ii2) = sig2(0,0) * Cjinv(ii1,ii2);
                    }
                }
                
        // # Term 1
        //  Eterm1 <- Eterm1 + wgtlev2[jj] * AfjArj[ 1:NX + (jj-1)*NX , ] %*% theta_rj.jj
        //  arma::mat Eterm1 = arma::zeros(NX,1);
        
        for (int ii=0;ii<NX;ii++){
        for (int kk=0;kk<NZ;kk++){ 
          Eterm1(ii,0) += wgtlev2[jj] * AfjArj( ii + jj*NX , kk ) * theta_rj_jj(kk,0);
                    }
                }
    
        // # Term 2
        // Eterm2 <- Eterm2 + wgtlev2[jj] * theta_rj.jj %*% t( theta_rj.jj) + sig2 * Cjinv
        // arma::mat Eterm2 = arma::zeros(NZ,NZ);    
        for (int ii1=0;ii1<NZ;ii1++){
        for (int ii2=0;ii2<NZ;ii2++){
            Eterm2(ii1,ii2) += wgtlev2[jj]* theta_rj_jj(ii1,0)*theta_rj_jj(ii2,0) + 
                           sig2(0,0)*Cjinv(ii1,ii2) ;
                       }
                 }
    
        //   # Term 3
        //   rj <- y[id.jj] - X[id.jj,] %*% theta  - Z[ id.jj , ] %*% theta_rj.jj
        //   A1 <- Cjinv %*% ArjArj[ 1:NZ + (jj-1)*NZ , ] 
        //   Eterm3 <- Eterm3 + sum( rj^2 * wgtlev1[ id.jj]  ) + sig2 * sum( diag( A1 ) )                 
        //   arma::mat Eterm3 = arma::zeros(1,1);
        //   A1(NZ,NZ)
    
        sumtr3 = 0 ;
        for (int ii=0;ii<NZ;ii++){
          A1(ii,ii) = 0;
          for (int kk=0;kk<NZ;kk++){
            A1(ii,ii) = A1(ii,ii) + Cjinv(ii,kk) * ArjArj( kk + jj*NZ , ii );
                    }
          sumtr3 += A1(ii,ii) * sig2(0,0) ; 
                }
        Eterm3(0,0) += sumtr3 ;   
    
        // compute residuals
        for (int nn=idcluster_table(jj,0);nn < idcluster_table(jj,1)+1;nn++){       
            rj = y[nn] - ypred(nn,0) ;
            for (int kk=0;kk<NZ;kk++){
               rj = rj - Z(nn,kk) * theta_rj_jj(kk,0) ;
                        }
            ypred_ranef(nn,0) = ypred(nn,0) + ( y[nn] - ypred(nn,0) ) - rj ;
            
            Eterm3 += wgtlev1[nn] * rj * rj ;               
                    }
                    
    // Rcpp::Rcout << "|" << std::endl ;              
        
        }  // end jj         
    //----

    //--- M step
    
    //    # theta
    //    theta <- XtX %*% ( Xty - res0$Eterm1 )
    arma::mat theta_new = arma::mat( XtX * Xty - XtX * Eterm1 ) ; 
    
    // # Tmat
    //    Tmat <- res0$Eterm2 / W2
    arma::mat Tmat_new = arma::zeros( NZ , NZ );
    double w2 = W2[0] ;
//    if ( reml == 1 ){
//        w2 = w2 - NZ ; 
//            }
    for (int ii=0;ii<NZ;ii++){
    for (int hh=0;hh<NZ;hh++){    
        Tmat_new(ii,hh) = Eterm2(ii,hh) / w2 ;
                }
            }    
    
    // # sig2
    //    sig2 <- ( res0$Eterm3 / W1 )[1,1]
    arma::mat sig2_new = arma::zeros(1,1);
    w2 = W1[0] ;
//    if ( reml == 1 ){
//        w2 = w2 - NX ;
//            }
    sig2_new(0,0) = Eterm3(0,0) / w2 ;
    
      return Rcpp::List::create( 
                _["theta"] = theta_new , _["Tmat"]=Tmat_new ,
                _["sig2"] = sig2_new ,
                _["theta_rj"] = theta_rj , _["Var_theta_rj"] = Var_theta_rj ,
                _["ypred_fixed"] = ypred , _["ypred_ranef"] = ypred_ranef
        ) ;  
}
//**********************************************************


//**********************************************************
// MLA2 sufficient statistics
Rcpp::List mla2_suffstat( Rcpp::NumericMatrix X , Rcpp::NumericMatrix Z , 
    Rcpp::NumericVector y , 
    Rcpp::NumericVector wgtlev2 , Rcpp::NumericVector wgtlev1 ,
    Rcpp::NumericVector wgttot , Rcpp::NumericMatrix idcluster_table ){

    int NZ = Z.ncol();
    int NX = X.ncol();
    int N = X.nrow();
    int NC = idcluster_table.nrow(); 

    //**************************************
    // XtX <- solve( t(X) %*% wgttotM %*% X )
    arma::mat Xt2=arma::zeros(NX,NX);  
    arma::mat Xty=arma::zeros(NX,1);
    
    for (int ii=0; ii<NX; ii++){
    for (int jj=0; jj<NX; jj++){
        for (int nn=0;nn<N;nn++){
           Xt2(ii,jj) += X(nn,ii)*X(nn,jj) * wgttot[nn] ;   
                    }
                }
            }
    arma::mat XtX = arma::inv(Xt2);

    // Xty <- t(X) %*% wgttotM %*% y
    for (int ii=0;ii<NX;ii++){
    for (int nn=0;nn<N;nn++){
        Xty(ii,0) += X(nn,ii) * y[nn] * wgttot[nn] ;
                }
            }
    
    // # A_rj^T A_rj
    // ArjArj <- matrix( 0 , nrow=NC*NZ , ncol=NZ)
    // AfjArj <- matrix( 0 , nrow=NC*NX , ncol=NZ)
    
    arma::mat ArjArj = arma::zeros( NC*NZ , NZ ); 
    arma::mat AfjArj = arma::zeros( NC*NX , NZ );
    
    // for (jj in 1:NC){
    //     id.jj <- seq( idcluster_table[jj,1] , idcluster_table[jj,2]  )
    //    ArjArj[ 1:NZ + (jj-1)*NZ , ] <- t( Z[ id.jj , ] ) %*% ( wgtlev1[ id.jj ] * Z[ id.jj , ]    )    
    //            }
    
    for (int jj=0;jj<NC; jj++){  // begin jj
    
    for (int ii1=0;ii1<NZ ;ii1++){
    for (int ii2=ii1; ii2<NZ  ; ii2++){
       for (int nn=idcluster_table(jj,0); nn < idcluster_table(jj,1)+1 ; nn++ ){
          ArjArj( ii1 + jj * NZ , ii2) += Z( nn , ii1 ) * Z(nn,ii2) * wgtlev1[ nn ] ;
                    }
            }
        }
    for (int ii2=0;ii2<NZ-1;ii2++){
    for (int ii1=ii2+1; ii1<NZ ; ii1++){
        ArjArj( ii1 + jj*NZ , ii2 ) = ArjArj( ii2 + jj*NZ , ii1 ) ;
                    }
                }

    // for (jj in 1:NC){
    //     id.jj <- seq( idcluster_table[jj,1] , idcluster_table[jj,2]  )   
    //    AfjArj[ 1:NX + (jj-1)*NX , ] <- t( X[ id.jj , ] ) %*% wgtlev1M[ id.jj , id.jj ] %*% Z[ id.jj , ]
    //            }
                
    for (int ii1=0;ii1<NX;ii1++){
    for (int ii2=0; ii2<NZ ; ii2++){
       for (int nn=idcluster_table(jj,0); nn < idcluster_table(jj,1)+1 ; nn++ ){
          AfjArj( ii1 + jj * NX , ii2) += X( nn , ii1 ) * Z(nn,ii2) * wgtlev1[ nn ] ;
                    }
            }
        }                                
        } // end jj
    
    // compute sum of weights        
    Rcpp::NumericVector W1(1);
    Rcpp::NumericVector W2(1);
    for (int jj=0;jj<NC;jj++){
       W2[0] += wgtlev2[jj] ;
            }
    for (int nn=0;nn<N;nn++){
       W1[0] += wgtlev1[nn] ;
            }        
                            
    return Rcpp::List::create( 
       _["ArjArj"] = ArjArj , _["AfjArj"] = AfjArj ,
       _["XtX"] = XtX , _["Xty"] = Xty ,
       _["W1"] = W1 , _["W2"] = W2 , _["Xt2"] = Xt2 
        ) ;  
    
        }
//**********************************************************
        
//**********************************************************
// MLA2 starting values
Rcpp::List mla2_inits( arma::mat Xa , Rcpp::NumericMatrix X ,
    Rcpp::NumericMatrix Z , 
    Rcpp::NumericVector y , int NZ , Rcpp::NumericVector wgtlev1 ,
    Rcpp::NumericVector wgttot ){
        
    // int NZ = Z.ncol();
    int NX = X.ncol();
    int N = X.nrow();
    // int NC = idcluster_table.nrow(); 

    //**************************************
    // XtX <- solve( t(X) %*% wgttotM %*% X )
    arma::mat Xt2=arma::zeros(NX,NX);  
    arma::mat Xty=arma::zeros(NX,1);
    
    for (int ii=0; ii<NX; ii++){
    for (int jj=0; jj<NX; jj++){
        for (int nn=0;nn<N;nn++){
           Xt2(ii,jj) += X(nn,ii)*X(nn,jj) * wgttot[nn] ;   
                    }
                }
            }
    arma::mat XtX = arma::inv(Xt2);
    
    // Xty <- t(X) %*% wgttotM %*% y
    for (int ii=0;ii<NX;ii++){
    for (int nn=0;nn<N;nn++){
        Xty(ii,0) += X(nn,ii) * y[nn] * wgttot[nn] ;
                }
            }        

    // compute sum of weights        
    Rcpp::NumericVector W1(1);
    for (int nn=0;nn<N;nn++){
       W1[0] += wgtlev1[nn] ;
            }                
            
    arma::mat theta = arma::mat( XtX * Xty ) ;    
    arma::mat ypred = arma::mat( Xa * theta ) ; 
    double var_res = 0 ;
    for (int nn=0;nn<N;nn++){
        var_res += pow( y[nn] - ypred(nn,0) , 2.0 ) * wgtlev1[nn] ;
               }
    var_res = var_res / W1[0] ;           
    
    arma::mat sig2 = arma::zeros(1,1);
    sig2(0,0) = 0.8 * var_res ;
    arma::mat Tmat = arma::zeros(NZ,NZ);
    Tmat(0,0) = 0.2*var_res ;
    if (NZ > 1 ){
        for (int vv=1;vv<NZ;vv++){
           Tmat(vv,vv) = 0.1*var_res ;
                    }
          }
      
    return Rcpp::List::create( 
       _["theta"] = theta , _["Tmat"] = Tmat , _["sig2"] =sig2       
        ) ;          
        }
//**********************************************************        


//**********************************************************
Rcpp::List bifie_mla2_estimation( arma::mat theta_init , arma::mat Tmat_init ,
        arma::mat sig2_init, int NX , int NZ , int NC , int N ,
        Rcpp::NumericMatrix X , Rcpp::NumericMatrix Z , 
        Rcpp::NumericVector y , Rcpp::NumericVector wgtlev2 , 
        Rcpp::NumericVector wgtlev1 ,
        Rcpp::NumericVector wgttot , Rcpp::NumericMatrix idcluster_table ,
        double globconv , int maxiter , Rcpp::NumericMatrix recov_constraint ,
        int is_rcov_constraint , int NRC ){
    
    
    int iter=0;
        
    //**********************************************
    // set inits
    arma::mat theta(NX,1);
    arma::mat Tmat(NZ,NZ);
    for (int ii=0;ii<NX;ii++){
        theta(ii,0) = theta_init(ii,0);
                }
    for (int ii=0;ii<NZ;ii++){
      for (int hh=0;hh<NZ;hh++){    
        Tmat(ii,hh) = Tmat_init(ii,hh);
                }
            }
    arma::mat sig2(1,1);
    sig2(0,0) = sig2_init(0,0);
    
    arma::mat Xa = rcppmat2armamat( X )["armamat"] ;    
    
  
        
    //***************************************
    // determine sufficient statistics ;        
    Rcpp::List res0 = mla2_suffstat(X ,  Z , y , wgtlev2 , wgtlev1 , wgttot , 
             idcluster_table );                
    arma::mat ArjArj = res0["ArjArj"];
    arma::mat AfjArj = res0["AfjArj"];
    arma::mat XtX = res0["XtX"];
    arma::mat Xty = res0["Xty"];
    arma::mat Xt2 = res0["Xt2"];
    Rcpp::NumericVector W1 = res0["W1"];
    Rcpp::NumericVector W2 = res0["W2"];           
    

    arma::mat sig20 = arma::zeros(1,1);
    arma::mat theta0 = arma::zeros(NX,1);
    arma::mat Tmat0 = arma::zeros(NZ,NZ);
    
//    Rcpp::List res22;
//    Rcpp::List res21;    

    double absval=1E4 ;
    
    Rcpp::List res2;
    
  //***************************** start algorithm 
   while ( ( absval > globconv ) & ( iter < maxiter ) ){    

    // define previous solutions
    sig20(0,0) = sig2(0,0);
    for (int ii=0;ii<NX;ii++){
        theta0(ii,0) = theta(ii,0) ;
                }
    for (int ii=0;ii<NZ;ii++){
      for (int hh=0;hh<NZ;hh++){    
        Tmat0(ii,hh) = Tmat(ii,hh);
                }
            }

    // E and M step
    res2 = mla2_emsteps( X, Z ,idcluster_table , Tmat , Xa , theta ,
        y , wgtlev2 , wgtlev1 , ArjArj , AfjArj , sig2 , W1 , W2 , XtX , Xty  ) ;
    arma::mat theta_new = res2["theta"];    
    for (int ii=0;ii<NX;ii++){
        theta(ii,0) = theta_new(ii,0) ;
                }
    arma::mat Tmat_new = res2["Tmat"];
    for (int ii=0;ii<NZ;ii++){
    for (int hh=0;hh<NZ;hh++){
        Tmat(ii,hh) = Tmat_new(ii,hh) ;
                }
            }
  //*** begin constraints random effects covariance matrix  
  // recov_constraint , is_rcov_constraint, NRC    
  if (is_rcov_constraint == 1){
      for (int hh=0;hh<NRC;hh++){
        Tmat(recov_constraint(hh,0), recov_constraint(hh,1))=recov_constraint(hh,2);
        Tmat(recov_constraint(hh,1), recov_constraint(hh,0))=recov_constraint(hh,2);      
      }
  }
  //*** end constraints    
    
  // include constraints for covariance matrix
  
    arma::mat sig2_new = res2["sig2"];
    sig2(0,0) = sig2_new(0,0);
        // Rcpp::Rcout << "sig2=" <<  sig2(0,0) <<  std::flush << std::endl ;
        // Rcpp::Rcout << "...................."  <<  std::flush << std::endl ;

    //--- check convergence    
    Rcpp::NumericVector absval2 = maxabsval_arma( theta , theta0 ) ;
    absval = absval2[0] ;
    // Rcpp::Rcout << "absval (theta)=" <<  absval <<  std::flush << std::endl ;
    absval2 = maxabsval_arma( Tmat , Tmat0 ) ;
    if (absval2[0] > absval ){
        absval = absval2[0] ;
                }    
    // Rcpp::Rcout << "absval (Tmat)=" <<  absval <<  std::flush << std::endl ;            
    absval2 = maxabsval_arma( sig2 , sig20 ) ;
    if (absval2[0] > absval ){
        absval = absval2[0] ;
                }
    // Rcpp::Rcout << "absval (sig2)=" <<  absval <<  std::flush << std::endl ;
                
    iter = iter + 1 ;
           // Rcpp::Rcout << "absval=" <<  absval <<  std::flush << std::endl ;    
    
        }
  // end EM algorithm
  //****************************************************
  
// Rcpp::Rcout << "iter" <<  iter << std::flush << std::endl ;  
  
    // extract posterior distribution random effects
    arma::mat theta_rj = res2["theta_rj"];
    arma::mat Var_theta_rj = res2["Var_theta_rj"];

    // post processing
    Rcpp::List res31 = mla2_postproc(  N , NX , NZ , y , X , Z , idcluster_table , 
           wgttot , theta , Tmat , sig2 , Var_theta_rj , AfjArj , wgtlev2 , Xt2) ;
    Rcpp::NumericVector vardec = res31["vardec"] ;
    arma::mat fvcov = res31["fvcov"] ;
    
    // collecting all parameters
    int NV = vardec.size() ; 
    int NP = NX + NZ*NZ + 1 + NV ; 
    
    // parameters to be saved.
    // theta   NX
    // Tmat    NZ x NZ => NZ*(NZ-1)/2
    // sig2    1
    // vardec  12 | variance decomposition        
    Rcpp::NumericVector pars_temp(NP);
    int vv = 0 ;
    // theta
    for (int ii=0;ii<NX;ii++){
        pars_temp[vv] = theta(ii,0) ;
        vv++;
        }
    // Tmat
    for (int ii=0;ii<NZ;ii++){
        for (int hh=ii;hh<NZ;hh++){
            pars_temp[vv] = Tmat(ii,hh) ;
            vv ++ ;
            if (ii<hh){
               pars_temp[vv] = Tmat( ii , hh ) / sqrt( Tmat(ii,ii) * Tmat(hh,hh) );
               vv ++ ;
                  }    
            }
        }                    
    // sig2
    pars_temp[vv] = sig2(0,0);
    vv ++ ;
    // vardec
    for (int ii=0;ii<NV;ii++){
         pars_temp[vv] = vardec[ii] ;
         vv ++ ;
                 }
    
    
    return Rcpp::List::create( 
        _["NX"] = NX , _["NZ"] = NZ , _["NC"] = NC , _["N"] = N ,
        _["W1"] = W1 , _["W2"] = W2 , _["fvcov"] = fvcov ,
        _["theta"] = theta , _["Tmat"] = Tmat , _["sig2"] = sig2 ,
        _["theta_rj"] = theta_rj  , _["Var_theta_rj"] = Var_theta_rj ,
        _["iter"] = iter  ,  _["vardec"] = vardec  
        ,    _["pars_"] = pars_temp   ,
            _["postproc"] = res31        ) ;
    }
//**********************************************************  

//**********************************************************
// estimation for replicate weights
Rcpp::List bifie_mla2_estimation_replicates( int N__ , int NC__ ,
    Rcpp::NumericVector wgttot__ , Rcpp::NumericMatrix wgtrep__ ,
    Rcpp::NumericVector wgtlev1__ , Rcpp::NumericVector wgtlev2__ ,
    Rcpp::NumericMatrix idcluster_table2 , arma::mat theta0 , arma::mat Tmat0 ,
    arma::mat sig20 , int NX , int NZ ,  Rcpp::NumericMatrix X__ ,  
    Rcpp::NumericMatrix Z__ , Rcpp::NumericVector y__ , 
    double globconv , int maxiter , int NP ,
    Rcpp::NumericMatrix recov_constraint , int is_rcov_constraint , int NRC
  ){

    Rcpp::NumericVector wgttot__rr(N__);
    Rcpp::NumericVector fac__rr(N__);
    Rcpp::NumericVector wgtlev2__rr(NC__);
    Rcpp::NumericVector wgtlev1__rr(N__);

    Rcpp::NumericVector wgtlev1a(N__);
    int RR = wgtrep__.ncol();
    double eps = 1E-8 ;
    double cl=0;
    int zz=0;
    
    
    Rcpp::NumericMatrix pars_temp(NP,RR);
    Rcpp::NumericVector iter_temp(RR);
    
    // begin rr
    for (int rr=0 ; rr <RR; rr++){

    
    for (int nn=0;nn<N__;nn++){
        fac__rr[nn] = wgtrep__(nn,rr) / wgttot__[nn] ;
        if (fac__rr[nn] < eps ){
            fac__rr[nn] = eps ;
                    }
        wgttot__rr[nn] = wgtrep__(nn,rr) ;
        if (wgttot__rr[nn] < eps ){
            wgttot__rr[nn] = eps ;
                    }

        wgtlev1__rr[nn] = wgtlev1__[nn] * fac__rr[nn] ;
                }
    
    cl=0;                
    for (int cc=0;cc<NC__;cc++){
         wgtlev2__rr[cc] = wgtlev2__[cc] * fac__rr[ idcluster_table2(cc,0) ] ;
         cl += wgtlev2__rr[cc] ;
                }
    for (int nn=0;nn<NC__;nn++){
          wgtlev2__rr[nn] = wgtlev2__rr[nn] * NC__ / cl ;
                }            
    for (int nn=0;nn<N__;nn++){
          wgttot__rr[nn] = wgttot__rr[nn] * NC__ / cl ;
                }            

    // rescaling weights
    wgtlev1a = rescale_lev1weights( idcluster_table2 , wgtlev1__rr );

    // estimation
    Rcpp::List res33 = bifie_mla2_estimation( theta0 , Tmat0 ,
         sig20 , NX , NZ ,  NC__ ,  N__ , X__ ,  Z__ , y__ , wgtlev2__rr , 
         wgtlev1a , wgttot__rr , idcluster_table2 , globconv , maxiter,
         recov_constraint, is_rcov_constraint, NRC
      ) ;
    Rcpp::NumericVector pars=res33["pars_"] ;
        for (int pp=0;pp<NP;pp++){
            pars_temp(pp,rr) = pars[pp] ;
                    }
    int iter3 = as<int>(res33["iter"]) ;
    iter_temp[rr] = iter3 ;
        
    zz ++ ;
    if ( zz == 9 ){
        zz = 0 ;
        Rcpp::Rcout << "-" <<    std::flush  ;    
        // Rcpp::Rcout << "-" <<  std::flush ;         
            }

    
    
    }  // end rr

    return Rcpp::List::create( 
        _["pars_temp"] = pars_temp ,
        _["iter_temp"] = iter_temp
            ) ;
}
//**********************************************************


