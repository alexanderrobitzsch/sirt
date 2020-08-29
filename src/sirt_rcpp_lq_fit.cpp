//// File Name: sirt_rcpp_lq_fit.cpp
//// File Version: 0.25



// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


// user includes


///********************************************************************
///** sirt_rcpp_lq_fit_analyze_matrix
// [[Rcpp::export]]
Rcpp::NumericMatrix sirt_rcpp_lq_fit_analyze_matrix( Rcpp::NumericMatrix X )
{
    int NR=X.nrow();
    int NC=X.ncol();
    int NZ=NR*NC;
    Rcpp::NumericMatrix Z(NZ,3);
    int zz=0;
    double eps=1e-15;
    for (int rr=0; rr<NR; rr++){
        for (int cc=0; cc<NC; cc++){
            if (std::abs(X(rr,cc))>eps){
                Z(zz,0)=rr;
                Z(zz,1)=cc;
                Z(zz,2)=X(rr,cc);
                zz++;
            }
        }
    }
    Z = Z( Rcpp::Range(0, zz-1), _);

    ///---- OUTPUT
    return Z;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_lq_fit_sparse_matrix_derivative
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_lq_fit_sparse_matrix_derivative(
    Rcpp::NumericMatrix Z, Rcpp::NumericVector h1, int px)
{
    int NZ=Z.nrow();
    Rcpp::NumericVector der(px);
    for (int zz=0; zz<NZ; zz++){
        der[ Z(zz,1) ] -= Z(zz,2)*h1[ Z(zz,0) ];
    }

    ///---- OUTPUT
    return der;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_lq_fit_matrix_mult
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_lq_fit_matrix_mult(
    Rcpp::NumericMatrix Z, Rcpp::NumericVector y, Rcpp::NumericVector beta)
{
    int N=y.size();
    int NZ=Z.nrow();
    Rcpp::NumericVector e(N);
    for (int nn=0; nn<N; nn++){
        e[nn] = y[nn];
    }
    for (int zz=0; zz<NZ; zz++){
        e[ Z(zz,0) ] -= Z(zz,2)*beta[ Z(zz,1) ];
    }

    ///---- OUTPUT
    return e;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_lq_fit_fct_optim
// [[Rcpp::export]]
double sirt_rcpp_lq_fit_fct_optim( Rcpp::NumericMatrix Z, Rcpp::NumericVector y,
    Rcpp::NumericVector beta, double pow, Rcpp::NumericVector w, double eps)
{
    Rcpp::NumericVector e=sirt_rcpp_lq_fit_matrix_mult(Z,y,beta);
    // val <- sum( w*(e^2 + eps )^(pow/2) )
    double val=0.0;
    int N=e.size();
    double pow2=pow/2.0;
    for (int nn=0; nn<N; nn++){
        val += w[nn]*std::exp( pow2*std::log( e[nn]*e[nn] + eps) );
    }

    ///---- OUTPUT
    return val;
}
///********************************************************************
