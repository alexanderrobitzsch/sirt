//// File Name: eigenvaluessirt.cpp
//// File Version: 4.04


// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


// [[packageincludes]]
#include "first_eigenvalue_sirt.h"
//#include "c:/Users/robitzsch/Dropbox/Eigene_Projekte/R-Routinen/IRT-Functions/sirt_package/1.15/sirt_work/src/first_eigenvalue_sirt__2.19.h"


///********************************************************************
///** eigenvaluesDsirt
// [[Rcpp::export]]
Rcpp::List eigenvaluesDsirt( Rcpp::NumericMatrix Xr, 
	int D , int maxit , double conv ){
 
     Rcpp::NumericVector d1(2) ;         
     double K=Xr.nrow() ;  
     Rcpp::NumericVector dvec(D) ;  
     arma::mat u(K,D) ;         
     Rcpp::List res2;  
       
     //**********  
     // set matrices  
     arma::mat X(Xr.begin(), K, K , false);   
       
     arma::mat X0 = X ;  
       
     for (int dd=0;dd<D;dd++){  
     // int dd = 0;  
     	// estimate first eigenvalue of reduced matrix  
     	res2 = firsteigenvalsirt(X0,maxit,conv,K);  
     	d1 = res2["lambda1"] ;  
      // Rcpp::Rcout << " d1 = " << d1[0] << std::endl ;	  
     	dvec[dd] = d1[0] ;  
     	arma::mat u1 = res2["u"] ;  
     	u.col(dd) = arma::mat( u1.col(0) ) ;  
          for (int ii1=0;ii1<K;ii1++){
          	X0(ii1,ii1) = X0(ii1,ii1) - dvec[dd] * u(ii1,dd)*u(ii1,dd) ;                     
            for (int ii2=ii1+1;ii2<K;ii2++){  
          	  X0(ii1,ii2) = X0(ii1,ii2) - dvec[dd] * u(ii1,dd)*u(ii2,dd) ;       
                X0(ii2,ii1) = X0(ii1,ii2) ; 
            }  
          }  
     	}  
                            
     ////////////////////////////////////  
     // OUTPUT:  
     return Rcpp::List::create(  
         Rcpp::_["d"]=dvec  ,   
         Rcpp::_["u"]= u   
            		) ;  
}




