//// File Name: bifiesurvey_rcpp_jack_data_prep.cpp
//// File Version: 1.26


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// #include <Rcpp.h>


using namespace Rcpp;
using namespace arma;



//*************************************************************************
//*************************************************************************
//*************************************************************************
//  generate Jackknife replicate weights TIMSS studies
//*************************************************************************
//*************************************************************************
//*************************************************************************

// ** bifiesurvey_rcpp_jackknife_timss
// [[Rcpp::export]]
Rcpp::NumericMatrix bifiesurvey_rcpp_jackknife_timss( Rcpp::NumericVector wgt,
        Rcpp::NumericVector jkzone, Rcpp::NumericVector jkrep,
        int RR, double jkfac, Rcpp::NumericVector prbar )
{
    int N = wgt.size();
    Rcpp::NumericMatrix wgtrep(N,RR);
    for (int rr=0;rr<RR;rr++){
        for (int nn=0;nn<N;nn++){
            if ( jkzone[nn] == rr ){
                wgtrep(nn,rr) = jkfac * wgt[nn] * jkrep[nn];
            } else {
                wgtrep(nn,rr) = wgt[nn];
            }
        }
        if (prbar[rr] == 1){
            Rcpp::Rcout << "-" <<  std::flush;
        }
    }
    return wgtrep;
}

    // Rcpp::Rcout << "tmp1 " <<  tmp1 <<  std::flush << std::endl;



//*************************************************************************
//*************************************************************************
//*************************************************************************
// generate bootstrap samples (IID assumption)
//*************************************************************************
//*************************************************************************
//*************************************************************************

//*** bifiesurvey_rcpp_bootstrap
// [[Rcpp::export]]
Rcpp::List bifiesurvey_rcpp_bootstrap( Rcpp::NumericVector cumwgt, Rcpp::NumericMatrix rand_wgt )
{
    int N = cumwgt.size();
    int WW = rand_wgt.ncol();
    arma::colvec wgt_temp(N);
    Rcpp::NumericMatrix wgtM(N,WW);
    int zz=0;
    int nn1=0;
    for (int ww=0;ww<WW;ww++){
        for (int nn=0;nn<N;nn++){
            wgt_temp(nn,0) = rand_wgt(nn,ww);
        }
        // sort entries
        arma::colvec wgt_tempsort = arma::sort( wgt_temp );
        zz=0;
        nn1=0;
        while (nn1<N){
            if (wgt_tempsort(nn1,0) < cumwgt[zz] ){
                wgtM(zz,ww) ++;
                nn1 ++;
            } else {
                zz ++;
            }
        }
    }

    //--- OUTPUT
    return Rcpp::List::create(
            Rcpp::Named("cumwgt") = cumwgt,
            Rcpp::Named("rand_wgt") = rand_wgt,
            Rcpp::Named("wgtM") = wgtM
        );
}

//*************************************************************************
//*************************************************************************
//*************************************************************************
// converting BIFIEdata object into BIFIEcdata object
//*************************************************************************
//*************************************************************************
//*************************************************************************

//-- bifie_bifiedata2bifiecdata
//++ bifiesurvey_rcpp_bifiedata2bifiecdata
// [[Rcpp::export]]
Rcpp::List bifiesurvey_rcpp_bifiedata2bifiecdata( Rcpp::NumericMatrix datalistM, int Nimp )
{
    int N = datalistM.nrow() / Nimp;
    int VV = datalistM.ncol();
    //************************************************
    //*** define indicator matrix
    Rcpp::NumericMatrix datalistM_ind(N, VV);
    double t1=0;
    // number of missing entries
    int Nmiss=0;
    for (int nn=0;nn<N;nn++){ // beg nn
        for (int vv=0;vv<VV;vv++){  // beg vv
            datalistM_ind(nn,vv) = 1;
            t1 = datalistM(nn,vv);
            if ( ! R_IsNA( datalistM(nn, vv ) ) ){
                for (int ii=1;ii<Nimp;ii++){ // beg ii
                    if ( datalistM( nn + ii*N, vv ) != t1 ){
                        datalistM_ind(nn,vv)=0;
                        Nmiss ++;
                        break;
                    }
                }  // end ii
            }  // end if R_isNA
        } // end vv
    } // end nn

    //************************************************
    // imputed data matrix
    // int ZZ=Nmiss*Nimp;
    // Rcpp::NumericMatrix datalistM_imputed(ZZ,4);
    // first column: imputed dataset
    // second column: subject case
    // third column: imputed variable
    // fourth column: value
    int ZZ=Nmiss;
    Rcpp::NumericMatrix datalistM_imputed(ZZ,Nimp);
    Rcpp::IntegerMatrix datalistM_impindex(ZZ,2);
    // first column: subject case
    // second column: imputed variable
    int zz=0;
    for (int vv=0;vv<VV;vv++){
        for (int nn=0;nn<N;nn++){
            if ( ! R_IsNA( datalistM(nn, vv ) ) ){
                if ( datalistM_ind(nn,vv) == 0 ){ // beg if missing entry
                    datalistM_impindex(zz,0)=nn;
                    datalistM_impindex(zz,1)=vv;
                    for (int ii=0;ii<Nimp;ii++){  // beg ii
                        //    datalistM_imputed(zz,0) = ii;
                        //    datalistM_imputed(zz,1) = nn;
                        //    datalistM_imputed(zz,2) = vv;
                        //    datalistM_imputed(zz,3) = datalistM( nn+ii*N, vv);
                        datalistM_imputed(zz,ii) = datalistM(nn+ii*N, vv );
                    } // end ii
                    zz ++;
                }    // end if missing entry
            }  // end R_IsNA
        } // end nn
    } // end vv

    //---- OUTPUT
    return Rcpp::List::create(
            Rcpp::Named("datalistM_ind") = datalistM_ind,
            Rcpp::Named("datalistM_imputed") = datalistM_imputed,
            Rcpp::Named("datalistM_impindex") = datalistM_impindex,
            Rcpp::Named("Nimp") = Nimp,
            Rcpp::Named("Nmiss") = Nmiss
        );
}



//*************************************************************************
//*************************************************************************
//*************************************************************************
// convert BIFIEcdata to BIFIEdata
//*************************************************************************
//*************************************************************************
//*************************************************************************

//-- bifie_bifiecdata2bifiedata
//** bifiesurvey_rcpp_bifiecdata2bifiedata
// [[Rcpp::export]]
Rcpp::List bifiesurvey_rcpp_bifiecdata2bifiedata( Rcpp::NumericMatrix datalistM_ind,
        Rcpp::NumericMatrix datalistM_imputed, int Nimp, Rcpp::NumericMatrix dat1,
        Rcpp::NumericMatrix datalistM_impindex )
{
    int N=dat1.nrow();
    int VV=dat1.ncol();
    int ZZ=Nimp*N;
    Rcpp::NumericMatrix datalistM(ZZ,VV);

    //**** non-imputed data
    for (int ii=0;ii<Nimp;ii++){   // beg ss
        for (int nn=0;nn<N;nn++){ // beg nn
            for (int vv=0;vv<VV;vv++){    // beg vv
                if ( datalistM_ind(nn,vv) == 1 ){ // beg non-imputed data
                    datalistM(nn+ii*N,vv) = dat1(nn,vv);
                }  // end non-imputed data
            } // end vv
        } // end nn
    } // end ss

    //**** imputed data
    //--- Rcpp::NumericMatrix datalistM_imputed(ZZ,4);
    //     first column: subject case
    //     second column: imputed variable
    int HH=datalistM_imputed.nrow();
    int nn_=0;
    int vv_=0;
    for ( int hh=0;hh<HH;hh++){
        nn_ = datalistM_impindex(hh,0);
        vv_ = datalistM_impindex(hh,1);
        for (int ii=0;ii<Nimp;ii++){
            datalistM( nn_ + ii*N, vv_ ) = datalistM_imputed(hh, ii);
        }
    }

    //---- OUTPUT
    return Rcpp::List::create(
            Rcpp::Named("datalistM") = datalistM,
            Rcpp::Named("Nimp") = Nimp
        );
}


//*************************************************************************
//*************************************************************************
//*************************************************************************
//*************************************************************************
// data preparation loading files with indicator matrix
//*************************************************************************
//*************************************************************************
//*************************************************************************

//-- bifie_bifiedata_stepwise
//++ bifiesurvey_rcpp_bifiedata_stepwise
// [[Rcpp::export]]
Rcpp::List bifiesurvey_rcpp_bifiedata_stepwise( Rcpp::NumericMatrix dat1,
        Rcpp::NumericMatrix dat_ind, int Nmiss )
{
    int N=dat1.nrow();
    int VV=dat1.ncol();
    //************************************************
    // imputed data matrix
    int ZZ=Nmiss;
    Rcpp::NumericMatrix datalistM_imputed(ZZ,4);
    // first column: imputed dataset
    // second column: subject case
    // third column: imputed variable
    // fourth column: value
    int zz=0;
    int ii=0;
    for (int vv=0;vv<VV;vv++){
        for (int nn=0;nn<N;nn++){
            if ( dat_ind(nn,vv) == 0 ){ // beg if missing entry
                datalistM_imputed(zz,0) = ii;
                datalistM_imputed(zz,1) = nn;
                datalistM_imputed(zz,2) = vv;
                datalistM_imputed(zz,3) = dat1( nn, vv);
                zz ++;
            }    // end if missing entry
        } // end nn
    } // end vv

    //--- OUTPUT
    return Rcpp::List::create(
            Rcpp::Named("datalistM_imputed") = datalistM_imputed
        );
}



