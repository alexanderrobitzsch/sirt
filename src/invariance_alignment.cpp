//// File Name: invariance_alignment.cpp
//// File Version: 1.02
//// File Last Change: 2017-02-17 10:52:09


// [[Rcpp::depends(RcppArmadillo)]]


// includes from the plugin
#include <RcppArmadillo.h>
#include <Rcpp.h>


using namespace Rcpp;


///********************************************************************
///** ia_optim_lambda
// [[Rcpp::export]]
Rcpp::NumericVector ia_optim_lambda( Rcpp::NumericMatrix lambda, 
	Rcpp::NumericVector psi0, Rcpp::NumericVector psi0b, 
	double align_scale, double align_pow, 
	Rcpp::NumericMatrix wgt, double eps, Rcpp::NumericMatrix group_combis ){
        
     int I = lambda.ncol() ;  
     int G = lambda.nrow() ;  
     Rcpp::NumericMatrix lambda1(G,I);  
     Rcpp::NumericMatrix lambda1b(G,I);       
     for (int ii=0;ii<I;ii++){  
     for (int gg=0;gg<G;gg++){  
     	lambda1(gg,ii) = lambda(gg,ii) / psi0[gg] ;  
             lambda1b(gg,ii) = lambda(gg,ii) / psi0b[gg] ;	  
     			}  
     		}  
       
     int GC = group_combis.nrow();  
     Rcpp::NumericVector fopt(GC);  
     Rcpp::NumericVector fopt1(GC);  
       
     for (int ii=0;ii<I;ii++){  
     for (int cc=0;cc<GC;cc++){  
     	//    fopt1 <- ( lambda1[ group.combis[,1] , ii ] - lambda1b[ group.combis[,2] , ii ] )^2  
     	fopt1[cc] = pow( lambda1( group_combis(cc,0) , ii) - lambda1b( group_combis(cc,1) , ii) , 2.0 ) ;  
     	//   fopt <- fopt + wgt[ group.combis[,1],ii] * wgt[ group.combis[,2],ii] *   
     	//                  ( fopt1 / align.scale^2 + eps )^align.pow  
     	fopt[cc] += wgt( group_combis(cc,0) , ii) * wgt( group_combis(cc,1) , ii) *  
     			   pow( fopt1[cc] / ( align_scale*align_scale ) + eps , align_pow ) ;  
     	}  
     }  
       
     // sum over the same indices  
     Rcpp::NumericVector res(G);  
     for (int cc=0;cc<GC;cc++){  
     	res[ group_combis(cc,0) ] += fopt[cc] ;		  
     	}  
     return res ;	  
}



///********************************************************************
///** ia_optim_nu
// [[Rcpp::export]]
Rcpp::NumericVector ia_optim_nu( Rcpp::NumericMatrix lambda, 
	Rcpp::NumericMatrix nu, Rcpp::NumericVector psi0_, 
	Rcpp::NumericVector psi0b, Rcpp::NumericVector alpha0, 
	Rcpp::NumericVector alpha0b, double align_scale , double align_pow , 
	Rcpp::NumericMatrix wgt, double eps , 
	Rcpp::NumericMatrix group_combis ){

     //    # optimization with respect to country SDs  
     //    nu1 <- nu - alpha0 * lambda   
     //    nu1b <- nu - alpha0b * lambda  
           
     int I = lambda.ncol() ;  
     int G = lambda.nrow() ;  
     Rcpp::NumericMatrix nu1(G,I);  
     Rcpp::NumericMatrix nu1b(G,I);  
       
     for (int ii=0;ii<I;ii++){  
	     for (int gg=0;gg<G;gg++){  
		nu1(gg,ii) = nu(gg,ii) - lambda(gg,ii) * alpha0[gg] ;  
		     nu1b(gg,ii) = nu(gg,ii) - lambda(gg,ii) * alpha0b[gg] ;	  
				}  
     		}  
       
     int GC = group_combis.nrow();  
     Rcpp::NumericVector fopt(GC);  
     Rcpp::NumericVector fopt1(GC);  
       
     for (int ii=0;ii<I;ii++){  
     for (int cc=0;cc<GC;cc++){  
     	//    fopt1 <- ( lambda1[ group.combis[,1] , ii ] - lambda1b[ group.combis[,2] , ii ] )^2  
     	fopt1[cc] = pow( nu1( group_combis(cc,0) , ii) - nu1b( group_combis(cc,1) , ii) , 2.0 ) ;  
     	//   fopt <- fopt + wgt[ group.combis[,1],ii] * wgt[ group.combis[,2],ii] *   
     	//                  ( fopt1 / align.scale^2 + eps )^align.pow  
     	fopt[cc] += wgt( group_combis(cc,0) , ii) * wgt( group_combis(cc,1) , ii) *  
     			   pow( fopt1[cc] / ( align_scale*align_scale ) + eps , align_pow ) ;  
     	}  
     }  
       
     // sum over the same indices  
     Rcpp::NumericVector res(G);  
     for (int cc=0;cc<GC;cc++){  
     	res[ group_combis(cc,0) ] += fopt[cc] ;		  
     }  
       
     return res ;	      
}

