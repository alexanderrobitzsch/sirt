//// File Name: sirt_rcpp_gom_em.cpp
//// File Version: 0.34




// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
// #include <Rcpp.h>

using namespace Rcpp;
using namespace arma;


///********************************************************************
///** sirt_rcpp_gom_em_calcpost
// [[Rcpp::export]]
Rcpp::List sirt_rcpp_gom_em_calcpost( Rcpp::IntegerMatrix DAT2,
    Rcpp::IntegerMatrix DAT2RESP, Rcpp::NumericMatrix PROBS,
    Rcpp::NumericMatrix DAT2IND, Rcpp::NumericVector PIK, Rcpp::NumericVector KK1,
    Rcpp::NumericVector weights)
{

    int N=DAT2.nrow();
    int I=DAT2.ncol();
    int TP=PROBS.ncol();
    int KK=KK1[0];

    //***  posterior distribution
    Rcpp::NumericVector pik1(TP);

    //---- calculate individual likelihood
    Rcpp::NumericMatrix fyiqk (N,TP);
    fyiqk.fill(1);
    for (int ii=0;ii<I;++ii){
        for (int nn=0;nn<N;++nn){
            if ( DAT2RESP(nn,ii)>0){
                for (int tt=0;tt<TP;++tt){
                    fyiqk(nn,tt) = fyiqk(nn,tt) * PROBS( 2*ii + DAT2(nn,ii), tt );
                }
            }
        }
    }

    //---- calculate posterior
    Rcpp::NumericMatrix fqkyi (N,TP);
    for (int nn=0;nn<N;++nn){
        double total = 0;
        for (int tt=0;tt<TP;++tt){
            fqkyi(nn,tt) = fyiqk(nn,tt)*PIK[tt];
            total += fqkyi(nn,tt);
        }
        for (int tt=0;tt<TP;++tt){
            fqkyi(nn,tt) = fqkyi(nn,tt)/total;
        }
    }
    //--- sum of weights
    double W = 0;
    for (int nn=0; nn<N; nn++){
        W += weights[nn];
    }

    //--- calculate counts
    for (int tt=0;tt<TP;++tt){
        pik1[tt] = 0;
        for (int nn=0;nn<N;++nn){
            pik1[tt] += weights[nn]*fqkyi(nn,tt);
        }
        pik1[tt] = pik1[tt] / W;
    }

    Rcpp::NumericMatrix nik (TP, I*(KK+1));
    Rcpp::NumericMatrix NIK (TP, I);
    for (int tt=0;tt<TP;++tt){
        for (int ii=0; ii < I; ++ii){
            for (int kk=0;kk<KK+1;++kk){
                for (int nn=0;nn<N;++nn){
                    nik( tt, ii + kk*I  ) += weights[nn]*DAT2IND(nn,ii+kk*I) *fqkyi(nn,tt);
                }  // end nn
                NIK(tt,ii) += nik(tt,ii+kk*I );
            } // end kk
        }  // end ii
    }  // end tt

    ///////////// O U T P U T   ///////////////////////////
    return Rcpp::List::create(
            Rcpp::Named("fyiqk") = fyiqk,
            Rcpp::Named("f.qk.yi") = fqkyi,
            Rcpp::Named("pi.k") = pik1,
            Rcpp::Named("n.ik") = nik,
            Rcpp::Named("N.ik") = NIK,
            Rcpp::Named("W") = W
        );
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_gom_em_likelihood
// [[Rcpp::export]]
Rcpp::NumericMatrix sirt_rcpp_gom_em_likelihood( Rcpp::NumericVector probs,
        int ncat, int TP, Rcpp::IntegerMatrix dat2, Rcpp::LogicalMatrix dat2resp)
{
    int N = dat2.nrow();
    int I = dat2.ncol();
    Rcpp::NumericMatrix fyiqk(N,TP);
    fyiqk.fill(1);
    for (int nn=0; nn<N; nn++){
        for (int ii=0; ii<I; ii++){
            if( dat2resp(nn,ii) ){
                for (int tt=0; tt<TP;tt++){
                    fyiqk(nn,tt) *= probs[ii+dat2(nn,ii)*I+tt*I*ncat];
                }
            }
        }
    }
    //--- output
    return fyiqk;
}
///********************************************************************


///********************************************************************
///** sirt_rcpp_gom_em_log_likelihood
// [[Rcpp::export]]
double sirt_rcpp_gom_em_log_likelihood( Rcpp::NumericMatrix fyiqk,
        Rcpp::NumericVector pik, Rcpp::NumericVector weights)
{
    int N = fyiqk.nrow();
    int TP = fyiqk.ncol();
    double ll = 0;
    double temp = 0;
    double eps = 1e-50;

    for (int nn=0; nn<N; nn++){
        temp = 0;
        for (int tt=0; tt<TP; tt++){
            temp += fyiqk(nn,tt)*pik[tt];
        }
        ll += weights[nn]*std::log(temp+eps);
    }
    ll = -2*ll;
    //--- output
    return ll;
}
///********************************************************************


///********************************************************************
///** sirt_rcpp_gom_em_loglike_gradient
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_gom_em_loglike_gradient( Rcpp::NumericVector probs,
        Rcpp::NumericVector probs_h, int ncat, int TP, Rcpp::IntegerMatrix dat2,
        Rcpp::LogicalMatrix dat2resp, Rcpp::NumericVector pik,
        Rcpp::IntegerVector items, Rcpp::NumericMatrix fyiqk,
        Rcpp::NumericVector weights)
{
    int IA = items.size();
    int I = dat2.ncol();
    int N = dat2.nrow();
    Rcpp::NumericVector ll_h(IA);
    int index_temp = 0;
    double ll_ind = 0;
    double fyiqk_temp = 0;
    double eps = 1e-50;

    for (int hh=0; hh<IA; hh++){
        int item_hh = items[hh];
        for (int nn=0; nn<N; nn++){
            ll_ind = 0;
            for (int tt=0; tt<TP; tt++){
                if (dat2resp(nn,item_hh)){
                    index_temp = item_hh + dat2(nn,item_hh)*I + tt*I*ncat;
                    fyiqk_temp = ( fyiqk(nn,tt) / probs[index_temp] ) * probs_h[index_temp];
                } else {
                    fyiqk_temp = fyiqk(nn,tt);
                }
                ll_ind += pik[tt]*fyiqk_temp;
            }
            ll_h[hh] += weights[nn]*std::log(ll_ind+eps);
        }
        ll_h[hh] = -2*ll_h[hh];
    }

    //--- output
    return ll_h;
}
///********************************************************************



