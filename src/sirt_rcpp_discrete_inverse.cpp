//// File Name: sirt_rcpp_discrete_inverse.cpp
//// File Version: 0.19


// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


// user includes

///********************************************************************
///** sirt_rcpp_discrete_inverse
// [[Rcpp::export]]
Rcpp::List sirt_rcpp_discrete_inverse( Rcpp::NumericVector x0,
        Rcpp::NumericVector y0, Rcpp::NumericVector y)
{
    int TP=x0.size();
    int N=y.size();
    Rcpp::IntegerVector ind(N);
    Rcpp::NumericVector x(N);
    int tt=0;
    bool iterate_nn;
    int ind_nn=0;
    int ind_nn1=0;
    double xdiff=x0[1]-x0[0];
    double xdiff2=xdiff/2.0;

    // start case nn
    for (int nn=0; nn<N; nn++){
        ind_nn=0;
        tt=0;
        iterate_nn=TRUE;
        while (iterate_nn){
            if (y[nn]>y0[tt]){
                ind_nn++;
                tt++;
            } else {
                iterate_nn=FALSE;
            }
            if (tt>TP-1){ iterate_nn=FALSE;}
        }
        ind[nn]=ind_nn;
        x[nn]=x0[ind_nn];
        if (ind_nn==0){
            ind[nn]=0;
            x[nn]=x0[ind_nn];
        } else {
            ind_nn1=ind_nn-1;
            x[nn]=x0[ind_nn]-xdiff2 + xdiff*(y[nn]-y0[ind_nn1])/(y0[ind_nn]-y0[ind_nn1]);
        }
    }

    //--- OUTPUT:
    return Rcpp::List::create(
            Rcpp::Named("x")=x,
            Rcpp::Named("ind")=ind
        );
}

