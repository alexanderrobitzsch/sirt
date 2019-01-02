//// File Name: sirt_rcpp_linking_haebara.cpp
//// File Version: 0.381


// [[Rcpp::depends(RcppArmadillo)]]


// includes from the plugin
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::interfaces(r, cpp)]]

using namespace Rcpp;



///********************************************************************
///** sirt_rcpp_linking_haebara_irf_2pl
// [[RcppNOexport]]
Rcpp::NumericVector sirt_rcpp_linking_haebara_irf_2pl( Rcpp::NumericVector theta,
    double a, double b)
{
    Rcpp::NumericVector x = a*(theta-b);
    Rcpp::NumericVector probs = Rcpp::plogis(x);
    //--- output
    return probs;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_linking_haebara_fct_optim_one_item
// [[Rcpp::export]]
double sirt_rcpp_linking_haebara_fct_optim_one_item( Rcpp::NumericVector theta,
    Rcpp::NumericVector prob_theta, Rcpp::NumericMatrix aM, Rcpp::NumericMatrix bM,
    Rcpp::NumericVector a, Rcpp::NumericVector b, Rcpp::NumericVector mu,
    Rcpp::NumericVector sigma, int ii, int ss, Rcpp::CharacterVector dist, double eps)
{
    double val=0;

    //** observed IRF
    Rcpp::NumericVector p_obs = sirt_rcpp_linking_haebara_irf_2pl( theta, aM(ii,ss), bM(ii,ss) );

    //** expected IRF
    //    a_exp <- a[ii] * sigma[ss]
    //    b_exp <- ( b[ii] - mu[ss] ) / sigma[ss]
    double a_exp = a[ii] * sigma[ss];
    double b_exp = ( b[ii] - mu[ss]) / sigma[ss];
    //    p_exp <- stats::plogis( a_exp * (theta - b_exp ) )
    Rcpp::NumericVector p_exp = sirt_rcpp_linking_haebara_irf_2pl( theta, a_exp, b_exp);

    //** compute distance
    Rcpp::NumericVector diff1 = p_obs - p_exp;
    diff1 = diff1*diff1;
    if (dist[0]=="L2"){
        val = Rcpp::sum( diff1*prob_theta );
    }
    if (dist[0]=="L1"){
        diff1 = Rcpp::sqrt( diff1 + eps);
        val = Rcpp::sum( diff1*prob_theta );
    }
    //--- output
    return val;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_linking_haebara_fct_optim
// [[Rcpp::export]]
double sirt_rcpp_linking_haebara_fct_optim( int NI, int NS, Rcpp::CharacterVector dist,
        Rcpp::NumericMatrix aM, Rcpp::NumericMatrix bM, Rcpp::NumericVector theta,
        Rcpp::NumericVector prob_theta, Rcpp::LogicalMatrix est_pars, Rcpp::NumericMatrix wgtM,
        Rcpp::NumericVector a, Rcpp::NumericVector b, Rcpp::NumericVector mu,
        Rcpp::NumericVector sigma, double eps )
{
    double val=0;
    double dist1;
    for (int ii=0; ii<NI; ii++){
        for (int ss=0; ss<NS; ss++){
            if (est_pars(ii,ss)){
                dist1 = sirt_rcpp_linking_haebara_fct_optim_one_item( theta,
                            prob_theta, aM, bM, a, b, mu, sigma, ii, ss, dist, eps );
                val += wgtM(ii,ss)*dist1;
            }  // end est_pars(ii,ss)
        }  // end ss
    } // end ii
    //--- output
    return val;
}
///********************************************************************


///********************************************************************
///** sirt_rcpp_linking_haebara_grad_optim_one_item
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_linking_haebara_grad_optim_one_item( Rcpp::NumericVector theta,
    Rcpp::NumericVector prob_theta, Rcpp::NumericMatrix aM, Rcpp::NumericMatrix bM,
    Rcpp::NumericVector a, Rcpp::NumericVector b, Rcpp::NumericVector mu,
    Rcpp::NumericVector sigma, int ii, int ss, Rcpp::CharacterVector dist, double eps,
    int NI, int NS, Rcpp::IntegerVector index_a, Rcpp::IntegerVector index_b,
    Rcpp::IntegerVector index_mu, Rcpp::IntegerVector index_sigma, Rcpp::NumericMatrix wgtM,
    Rcpp::NumericVector grad0)
{
    Rcpp::NumericVector grad = Rcpp::clone(grad0);
    Rcpp::NumericVector der_basis, der_t, diff2, dist2;
    int ind;
    double der_t1;

    //* observed IRF
    Rcpp::NumericVector p_obs = sirt_rcpp_linking_haebara_irf_2pl( theta, aM(ii,ss), bM(ii,ss) );
    //* expected IRF
    double a_exp = a[ii] * sigma[ss];
    double b_exp = ( b[ii] - mu[ss] ) / sigma[ss];
    Rcpp::NumericVector p_exp = sirt_rcpp_linking_haebara_irf_2pl( theta, a_exp, b_exp );
    Rcpp::NumericVector der = p_exp*(1-p_exp);
    diff2 = p_obs - p_exp;
    if (dist[0]=="L2"){
        der_basis = -2*diff2*prob_theta*der;
    }
    if (dist[0]=="L1"){
        dist2 = diff2*diff2;
        der_basis = - diff2*prob_theta*der / Rcpp::sqrt(dist2+eps);
    }
    double w_t = wgtM(ii,ss);
    //- a
    der_t = sigma[ss]*theta + mu[ss] - b[ii];
    ind = index_a[ii];
    grad[ind] += w_t*Rcpp::sum(der_basis*der_t);
    //- b
    der_t1 = -a[ii];
    ind = index_b[ii];
    grad[ind] += w_t*Rcpp::sum(der_basis*der_t1 );
    //- mu
    if (ss>0){
        der_t1 = a[ii];
        ind = index_mu[ss-1];
        grad[ind] += w_t*Rcpp::sum(der_basis*der_t1 );
    }
    //- sigma
    if (ss>0){
        der_t = a[ii]*theta;
        ind = index_sigma[ss-1];
        grad[ind] += w_t*Rcpp::sum(der_basis*der_t );
    }

    //--- output
    return grad;
}
///********************************************************************



///********************************************************************
///** sirt_rcpp_linking_haebara_grad_optim
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_linking_haebara_grad_optim( int NI, int NS,
    Rcpp::CharacterVector dist, Rcpp::NumericMatrix aM, Rcpp::NumericMatrix bM,
    Rcpp::NumericVector theta, Rcpp::NumericVector prob_theta, Rcpp::LogicalMatrix est_pars,
    Rcpp::NumericMatrix wgtM, Rcpp::NumericVector a, Rcpp::NumericVector b, Rcpp::NumericVector mu,
    Rcpp::NumericVector sigma, double eps, Rcpp::IntegerVector index_a, Rcpp::IntegerVector index_b,
    Rcpp::IntegerVector index_mu, Rcpp::IntegerVector index_sigma, Rcpp::CharacterVector parnames,
    int NP )
{
    Rcpp::NumericVector grad(NP, 0.0);
    for (int ii=0; ii<NI; ii++){
        for (int ss=0; ss<NS; ss++){
            if (est_pars(ii,ss)){
                grad = sirt_rcpp_linking_haebara_grad_optim_one_item( theta,
                            prob_theta, aM,  bM, a, b, mu, sigma, ii, ss, dist, eps,
                            NI, NS, index_a, index_b, index_mu, index_sigma,
                            wgtM, grad);
            }  // end est_pars(ii,ss)
        }  // end ss
    } // end ii

    //--- output
    return grad;
}
///********************************************************************

