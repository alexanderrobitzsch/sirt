//// File Name: sirt_rcpp_xxirt.cpp
//// File Version: 0.363



// [[Rcpp::depends(RcppArmadillo)]]

// #include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


///********************************************************************
///** sirt_rcpp_xxirt_compute_posterior_expected_counts
// [[Rcpp::export]]
Rcpp::NumericMatrix sirt_rcpp_xxirt_compute_posterior_expected_counts(
        Rcpp::LogicalMatrix dat1_resp_gg, Rcpp::NumericMatrix p_aj_xi_gg,
        Rcpp::NumericVector weights_gg)
{
    int N = dat1_resp_gg.nrow();
    int I = dat1_resp_gg.ncol();
    int TP = p_aj_xi_gg.ncol();

    Rcpp::NumericMatrix nij(I, TP);
    double val=0;

    //*** loop over items and categories
    for (int ii=0;ii<I;ii++){
        for (int tt=0;tt<TP;tt++){
            val=0;
            for (int nn=0;nn<N;nn++){
                if ( dat1_resp_gg(nn,ii) ){
                    val +=  weights_gg(nn)*p_aj_xi_gg(nn,tt);
                }
            }  // end nn
            nij(ii,tt) = val;
        }   // end tt
    }   // end ii

    //--- OUTPUT
    return nij;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_xxirt_compute_likelihood
// [[Rcpp::export]]
Rcpp::NumericMatrix sirt_rcpp_xxirt_compute_likelihood(
        Rcpp::IntegerMatrix dat, Rcpp::LogicalMatrix dat_resp_bool,
        Rcpp::NumericMatrix probs, int TP, int maxK )
{
    int N = dat.nrow();
    int I = dat.ncol();
    Rcpp::NumericMatrix p_xi_aj(N, TP);

    for (int nn=0;nn<N;nn++){
        for (int tt=0;tt<TP;tt++){
            p_xi_aj(nn,tt) = 1;
        }
        for (int ii=0;ii<I;ii++){
            if ( dat_resp_bool(nn,ii) ){
                for (int tt=0;tt<TP;tt++){
                    p_xi_aj(nn,tt) = p_xi_aj(nn,tt) * probs(ii, dat(nn,ii) + tt*maxK );
                }
            }
        }
    }

    //---- OUTPUT
    return p_xi_aj;
}
///********************************************************************



///********************************************************************
///** sirt_rcpp_xxirt_hessian_reduced_probs
// [[Rcpp::export]]
Rcpp::NumericMatrix sirt_rcpp_xxirt_hessian_reduced_probs(
        Rcpp::IntegerMatrix dat, Rcpp::LogicalMatrix dat_resp_bool,
        Rcpp::NumericMatrix probs_ratio, int TP, int maxK,
        int itemnr, int itemnr2, bool use_itemnr2, Rcpp::NumericMatrix p_xi_aj )
{
    int N = dat.nrow();
    Rcpp::NumericMatrix t1(N, TP);

    int ii=itemnr;

    for (int nn=0;nn<N;nn++){
        if (dat_resp_bool(nn,itemnr)){
            for (int tt=0; tt<TP; tt++){
                t1(nn,tt) = probs_ratio(ii, dat(nn,ii) + tt*maxK );
            }
        } else {
            for (int tt=0; tt<TP; tt++){
                t1(nn,tt) = 1;
            }
        }
    }
    if (use_itemnr2){
        ii = itemnr2;
        for (int nn=0;nn<N;nn++){
            if (dat_resp_bool(nn,itemnr2)){
                for (int tt=0; tt<TP; tt++){
                    t1(nn,tt) = t1(nn,tt) * probs_ratio(ii, dat(nn,ii) + tt*maxK );
                }
            }
        }
    }

    for (int nn=0;nn<N;nn++){
        for (int tt=0; tt<TP; tt++){
            t1(nn,tt) = t1(nn,tt)*p_xi_aj(nn,tt);
        }
    }

    //---- OUTPUT
    return t1;
}
///********************************************************************



///********************************************************************
///** sirt_rcpp_xxirt_newton_raphson_reduced_probs
// [[Rcpp::export]]
double sirt_rcpp_xxirt_newton_raphson_derivative_par(
        Rcpp::IntegerMatrix dat, Rcpp::LogicalMatrix dat_resp_bool,
        Rcpp::NumericMatrix ratio, Rcpp::NumericMatrix p_xi_aj,
        int item, Rcpp::NumericMatrix prior_Theta,
        Rcpp::IntegerVector group0, Rcpp::NumericVector weights,
        Rcpp::NumericVector ll_case0, double eps)
{
    int N = dat.nrow();
    int TP = p_xi_aj.ncol();
    double temp = 0;
    int ii = item-1;
    double ll_case_der_nn=0;
    double grad_pp=0;

    for (int nn=0; nn<N; nn++){
        if (dat_resp_bool(nn,ii)){
            ll_case_der_nn=0;
            for (int tt=0; tt<TP; tt++){
                temp=p_xi_aj(nn,tt)*ratio( dat(nn,ii), tt);
                ll_case_der_nn += prior_Theta(tt,group0[nn])*temp;
            }
            grad_pp -= weights[nn] * ll_case_der_nn / ( ll_case0[nn] + eps );
        }
    }

    //---- OUTPUT
    return grad_pp;
}
///********************************************************************


// if ( ! R_IsNA( resp(nn,ii) ) ){
