//// File Name: sirt_rcpp_inference_jackknife.h
//// File Version: 0.01

#ifndef _SIRT_SIRT_RCPP_INFERENCE_JACKKNIFE_H
#define _SIRT_SIRT_RCPP_INFERENCE_JACKKNIFE_H
 
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

Rcpp::List sirt_rcpp_inference_jackknife( Rcpp::NumericMatrix PARS );

#endif // _SIRT_SIRT_RCPP_INFERENCE_JACKKNIFE_H
