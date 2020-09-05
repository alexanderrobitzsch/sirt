//// File Name: sirt_rcpp_ccov.cpp
//// File Version: 0.17



// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


// user includes


///********************************************************************
///** sirt_rcpp_ccov_np_compute_ccov_sum_score
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_ccov_np_compute_ccov_sum_score(
        Rcpp::IntegerVector index, int NS, Rcpp::NumericMatrix data )
{
    Rcpp::NumericVector ccov(NS);
    Rcpp::NumericVector m1(NS);
    Rcpp::NumericVector m2(NS);
    Rcpp::NumericVector count(NS);
    int N=data.nrow();
    double eps=1e-15;
    int index_nn=0;
    for (int nn=0; nn<N; nn++){
        index_nn=index[nn];
        count[index_nn] ++;
        m1[index_nn] += data(nn,0);
        m2[index_nn] += data(nn,1);
        ccov[index_nn] += data(nn,0)*data(nn,1);
    }
    for (int ss=0; ss<NS; ss++){
        count[ss] += eps;
        ccov[ss] = ccov[ss] / count[ss] - m1[ss] * m2[ss] / ( count[ss] * count[ss] );
        //* unbiased estimation of conditional covariance
        // ccov[ss] = ( ccov[ss] - m1[ss] * m2[ss] / count[ss] ) / ( count[ss] - 1 );
    }

    ///---- OUTPUT
    return ccov;
}
///********************************************************************
