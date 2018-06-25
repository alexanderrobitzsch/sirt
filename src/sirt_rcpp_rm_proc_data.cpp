//// File Name: sirt_rcpp_rm_proc_data.cpp
//// File Version: 0.25


// [[Rcpp::depends(RcppArmadillo)]]

// #include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
// using namespace arma;


///********************************************************************
///** sirt_rcpp_rm_proc_datasets_na_indicators
// [[Rcpp::export]]
Rcpp::List sirt_rcpp_rm_proc_datasets_na_indicators(Rcpp::NumericMatrix dat, int K)
{
    int N = dat.nrow();
    int I = dat.ncol();
    Rcpp::IntegerMatrix dat_resp(N,I);
    Rcpp::NumericMatrix dat2(N,I);
    dat_resp.fill(0);
    dat2.fill(0);

    int K1 = K + 1;
    int NR = N * I * K1;
    int NI = N * I;
    Rcpp::NumericVector dat2_ind_resp(NR);

    for (int nn=0; nn<N; nn++){
        for (int ii=0; ii<I; ii++){
            if ( ! R_IsNA( dat(nn,ii) ) ){
                dat_resp(nn,ii) = 1;
                dat2(nn,ii) = dat(nn,ii);
                dat2_ind_resp[ nn + ii*N + dat(nn,ii)*NI ] = 1;
            }
        }
    }

    //-- output
    return Rcpp::List::create(
                Rcpp::Named("dat2") = dat2,
                Rcpp::Named("dat_resp") = dat_resp,
                Rcpp::Named("dat2_ind_resp") = dat2_ind_resp
            );
}
///********************************************************************


///********************************************************************
///** sirt_rcpp_rm_proc_expand_dataset
// [[Rcpp::export]]
Rcpp::NumericMatrix sirt_rcpp_rm_proc_expand_dataset(Rcpp::NumericMatrix dat,
        Rcpp::IntegerVector rater0,    Rcpp::IntegerVector pid0, int N, int R)
{
    int ND = dat.nrow();
    int I = dat.ncol();
    Rcpp::NumericMatrix dat2(N,I*R);
    dat2.fill(NA_REAL);
    int rater_dd=0;
    int IR = 0;

    for (int dd = 0; dd<ND; dd++){
        rater_dd = rater0[dd];
        IR = I * rater_dd;
        for (int ii = 0; ii<I; ii++){
            dat2( pid0[dd], IR + ii) = dat(dd,ii);
        }
    }

    //-- output
    return dat2;
}
///********************************************************************
