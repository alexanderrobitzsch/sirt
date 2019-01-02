//// File Name: sirt_rcpp_eigenvalues.cpp
//// File Version: 4.28


// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;



///********************************************************************
///** sirt_rcpp_first_eigenvalue
// [[Rcpp::export]]
Rcpp::List sirt_rcpp_first_eigenvalue(arma::mat X, int maxit, double conv, double K)
{
    double lambda=0;
    double lambdadiff=1000;
    double lambda_old, lambda_temp;
    Rcpp::NumericVector lambda_est(2);
    int iter = 0;

    //**********
    // set matrices
    arma::mat Xz;
    arma::colvec z(K);
    double temp1 = 1 / std::sqrt(K);
    for (int ii=0; ii<K; ii++){
        z[ii] = temp1;
    }

    ///---- algorithm
    while ( ( iter < maxit ) & ( lambdadiff > conv) ){
        lambda_old = lambda;
        Xz = arma::mat( X * z );
        lambda_temp = 0;
        for (int ii=0; ii<K; ii++){
            lambda_temp += Xz[ii]*Xz[ii];
        }
        lambda = std::sqrt(lambda_temp);
        lambdadiff = std::abs(lambda - lambda_old);
        z = Xz / lambda;
        iter ++;
    }
    lambda_est[0] = lambda;

    //---- OUTPUT:
    return Rcpp::List::create(
                Rcpp::Named("u")=z,
                Rcpp::Named("lambda1")=lambda_est
            );
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_D_eigenvalues
// [[Rcpp::export]]
Rcpp::List sirt_rcpp_D_eigenvalues( Rcpp::NumericMatrix Xr, int D, int maxit, double conv )
{
    Rcpp::NumericVector d1(2);
    double K=Xr.nrow();
    Rcpp::NumericVector dvec(D);
    arma::mat u(K,D);
    Rcpp::List res2;

    //******* set matrices
    arma::mat X0(Xr.begin(), K, K, false);
    for (int dd=0; dd<D; dd++){
        // estimate first eigenvalue of reduced matrix
        res2 = sirt_rcpp_first_eigenvalue(X0, maxit, conv, K);
        d1 = res2["lambda1"];
        dvec[dd] = d1[0];
        arma::mat u1 = res2["u"];
        u.col(dd) = arma::mat( u1.col(0) );
        for (int ii1=0; ii1<K; ii1++){
            X0(ii1,ii1) = X0(ii1,ii1) - dvec[dd] * u(ii1,dd)*u(ii1,dd);
            for (int ii2=ii1+1; ii2<K; ii2++){
                X0(ii1,ii2) = X0(ii1,ii2) - dvec[dd] * u(ii1,dd)*u(ii2,dd);
                X0(ii2,ii1) = X0(ii1,ii2);
            }
        }
    }

    //----- OUTPUT:
    return Rcpp::List::create(
                Rcpp::Named("d")=dvec,
                Rcpp::Named("u")=u
            );
}
///********************************************************************

// Rcpp::Rcout << " d1 = " << d1[0] << std::endl;

