//// File Name: sirt_rcpp_pbvnorm.h
//// File Version: 0.03

#ifndef _SIRT_SIRT_RCPP_PBVNORM_H
#define _SIRT_SIRT_RCPP_PBVNORM_H
 
#include <Rcpp.h>
using namespace Rcpp;

double sirt_rcpp_pnorm0( double z);

Rcpp::NumericVector sirt_rcpp_pnorm( Rcpp::NumericVector x);

double sirt_rcpp_pbvnorm0( double h1, double hk, double r);

Rcpp::NumericVector sirt_rcpp_pbvnorm( Rcpp::NumericVector x, Rcpp::NumericVector y,
        Rcpp::NumericVector rho);

double sirt_rcpp_dbvnorm0( double x, double y, double rho, bool use_log);

Rcpp::NumericVector sirt_rcpp_dbvnorm( Rcpp::NumericVector x, Rcpp::NumericVector y,
        Rcpp::NumericVector rho, bool use_log);

#endif // _SIRT_SIRT_RCPP_PBVNORM_H
