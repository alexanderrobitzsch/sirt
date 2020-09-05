//// File Name: sirt_rcpp_pbvnorm.cpp
//// File Version: 0.05


// [[Rcpp::depends(RcppArmadillo)]]

// #include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[RcppNOinterfaces(r, cpp)]]


const double pi = 3.1415926535897;

///********************************************************************
///** sirt_rcpp_pnorm0
// [[Rcpp::export]]
double sirt_rcpp_pnorm0( double z)
{
    double y = ::Rf_pnorm5(z, 0.0, 1.0, 1, 0);
    //--- OUTPUT
    return y;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_pnorm
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_pnorm( Rcpp::NumericVector x)
{
    int N = x.size();
    Rcpp::NumericVector y(N);
    for (int nn=0; nn<N; nn++){
        y[nn] = sirt_rcpp_pnorm0( x[nn] );
    }
    //--- OUTPUT
    return y;
}
///********************************************************************

///********************************************************************
///** Drezner & Wesolowksy, 1990, JCSC
///** sirt_rcpp_pbvnorm0
// [[Rcpp::export]]
double sirt_rcpp_pbvnorm0( double h1, double hk, double r)
{
    int NX=5;
    Rcpp::NumericVector X(NX);
    Rcpp::NumericVector W(NX);
    // data
    X[0]=.04691008;
    X[1]=.23076534;
    X[2]=.5;
    X[3]=.76923466;
    X[4]=.95308992;
    W[0]=.018854042;
    W[1]=.038088059;
    W[2]=.0452707394;
    W[3]=.038088059;
    W[4]=.018854042;
    // declarations
    double bv = 0;
    double r1, r2, rr, rr2, r3, h3, h5, h6, h7, aa, ab, h11;
    double cor_max = 0.7;
    double bv_fac1 = 0.13298076;
    double bv_fac2 = 0.053051647;

    // computation
    double h2 = hk;
    double h12 = (h1*h1+h2*h2)/2;
    double r_abs = std::abs(r);
    if (r_abs > cor_max){
        r2 = 1.0 - r*r;
        r3 = std::sqrt(r2);
        if (r<0){
            h2 = -h2;
        }
        h3 = h1*h2;
        h7 = std::exp( -h3 / 2.0);
        if ( r_abs < 1){
            h6 = std::abs(h1-h2);
            h5 = h6*h6 / 2.0;
            h6 = h6 / r3;
            aa = 0.5 - h3 / 8.0;
            ab = 3.0 - 2.0 * aa * h5;
            bv = bv_fac1*h6*ab*(1-sirt_rcpp_pnorm0(h6))-std::exp(-h5/r2)*(ab + aa*r2)*bv_fac2;
            for (int ii=0; ii<NX; ii++){
                r1 = r3*X[ii];
                rr = r1*r1;
                r2 = std::sqrt( 1.0 - rr);
                bv += - W[ii]*std::exp(- h5/rr)*(std::exp(-h3/(1.0+r2))/r2/h7 - 1.0 - aa*rr);
            }
        }
        h11 = std::min(h1,h2);
        bv = bv*r3*h7 + sirt_rcpp_pnorm0(h11);
        if (r < 0){
            bv = sirt_rcpp_pnorm0(h1) - bv;
        }
    } else {
        h3=h1*h2;
        for (int ii=0; ii<NX; ii++){
            r1 = r*X[ii];
            rr2 = 1.0 - r1*r1;
            bv += W[ii] * std::exp(( r1*h3 - h12)/rr2)/ std::sqrt(rr2);
        }
        bv = sirt_rcpp_pnorm0(h1)*sirt_rcpp_pnorm0(h2) + r*bv;
    }
    //--- OUTPUT
    return bv;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_pbvnorm
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_pbvnorm( Rcpp::NumericVector x, Rcpp::NumericVector y,
            Rcpp::NumericVector rho)
{
    int N = x.size();
    Rcpp::NumericVector res(N);
    for (int ii=0; ii<N; ii++){
        res[ii] = sirt_rcpp_pbvnorm0(x[ii], y[ii], rho[ii]);
    }
    //--- OUTPUT
    return res;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_dbvnorm0
// [[Rcpp::export]]
double sirt_rcpp_dbvnorm0( double x, double y, double rho, bool use_log)
{
    double pi2 = 2*pi;
    double r2 = 1-rho*rho;
    double r3 = std::sqrt(r2);
    double z = x*x - 2*rho*x*y + y*y;
    z = - z / r2 / 2.0;
    if ( ! use_log ){
        z = std::exp(z) / pi2 / r3;
    } else {
        z += - std::log(r3*pi2);
    }
    //--- OUTPUT
    return z;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_dbvnorm
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_dbvnorm( Rcpp::NumericVector x, Rcpp::NumericVector y,
            Rcpp::NumericVector rho, bool use_log)
{
    int N = x.size();
    Rcpp::NumericVector res(N);
    for (int ii=0; ii<N; ii++){
        res[ii] = sirt_rcpp_dbvnorm0(x[ii], y[ii], rho[ii], use_log);
    }
    //--- OUTPUT
    return res;
}
///********************************************************************
