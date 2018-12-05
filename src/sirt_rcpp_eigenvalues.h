//// File Name: sirt_rcpp_eigenvalues.h
//// File Version: 4.25

#ifndef _SIRT_SIRT_RCPP_EIGENVALUES_H
#define _SIRT_SIRT_RCPP_EIGENVALUES_H
 
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

Rcpp::List sirt_rcpp_first_eigenvalue(arma::mat X, int maxit, double conv, double K);

Rcpp::List sirt_rcpp_D_eigenvalues( Rcpp::NumericMatrix Xr, int D, int maxit, double conv );

#endif // _SIRT_SIRT_RCPP_EIGENVALUES_H
