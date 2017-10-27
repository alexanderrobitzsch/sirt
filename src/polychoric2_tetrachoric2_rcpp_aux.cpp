//// File Name: polychoric2_tetrachoric2_rcpp_aux.cpp
//// File Version: 3.07


// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


// [[packageincludes]]
#include "pbivnorm_rcpp_aux.h"
// #include "c:/Users/robitzsch/Dropbox/Eigene_Projekte/R-Routinen/IRT-Functions/sirt_package/1.15/sirt_work/src/pbivnorm_rcpp_aux__0.52.h"



///********************************************************************
///** tetrachoric2_rcpp_aux
// [[Rcpp::export]]
Rcpp::NumericVector tetrachoric2_rcpp_aux( Rcpp::NumericMatrix dfr, 
	double numdiffparm, int maxiter ){

     int ZZ = dfr.nrow() ;   
       
     Rcpp::NumericVector x(1);  
     Rcpp::NumericVector y(1);  
     Rcpp::NumericVector rho(1);  
     Rcpp::NumericVector p11(1);  
     Rcpp::NumericVector rhores(ZZ);  
       
     Rcpp::NumericVector L0(1) , L0h(1) ;  
     double L1;  
     double aincr=1 ;  
     double incr=0;  
     double maxincr = 1e-9 ;
     double eps=1e-10;
     
     int iter=0 ;  
       
     for (int zz = 0 ; zz<ZZ;zz++){  
     	  
     	iter=0 ;   
     	aincr=1 ;  
     	  
     	x[0] = dfr(zz,10 );  
     	y[0] = dfr(zz,11 );  
     	rho[0] = dfr(zz,15) ;  
     	p11[0] = dfr(zz,7) ;  
     	  
     	while ( ( iter < maxiter ) & ( aincr > maxincr ) ){  
     		L0 = pbivnorm2_rcpp( x , y , rho ) ;  
     		L0h = pbivnorm2_rcpp( x , y , rho+numdiffparm ) ;  
     		L1 = ( L0h[0] - L0[0] ) / numdiffparm ;  
     		// Newton-Raphson step according Divgi (1979)  
     		incr = ( L0[0] - p11[0] ) / ( L1 + eps ) ;  
     		rho[0] = rho[0] - incr ;  
     		iter ++ ;  
     		if (incr <0){ aincr = - incr ; } else { aincr=incr ; }  
     //		Rcout << "*** zz " << zz << " iter " << iter << "  incr " << incr << std::endl;  
     //		Rcout << "rho " << rho[0] << std::endl;  
     			}  
     	rhores[zz] = rho[0] ;  
     }  
       
     return rhores ;  
}


     // Rcout << "prob1 " << prob1 << std::endl;

     
///********************************************************************
///** polychoric2_aux_rcpp
// [[Rcpp::export]]
Rcpp::List polychoric2_aux_rcpp( Rcpp::NumericMatrix dat, 
	int maxK , int maxiter ){
      
     // int N = dat.nrow() ;  
     int I = dat.ncol() ;  
     int maxK3 = maxK + 3 ;  
       
     Rcpp::NumericVector v1 , tmp1 , tmp2 ;  
     Rcpp::NumericVector v2;  
     Rcpp::List res ;  
       
     Rcpp::NumericMatrix rho(I,I);  
     Rcpp::NumericMatrix Nobs(I,I);  
     Rcpp::NumericMatrix thresh(I,maxK3) ;  
     Rcpp::NumericVector maxcat(I);  
     Rcpp::NumericVector tmp3;  
     Rcpp::NumericVector Ntot_used(I);  
       
     for (int ii=0 ; ii < I-1;ii++){  
     for (int jj=ii+1;jj<I;jj++){  
     	v1 = dat(_,ii) ;  
     	v2 = dat(_,jj) ;	  
     	res=polychoric2_itempair( v1 , v2 , maxK , maxiter ) ;  
     	tmp1 = res["rho"] ;  
     	rho(ii,jj) = tmp1[0] ;  
     	rho(jj,ii) = rho(ii,jj);  
     	tmp1 = res["Ntotal"] ;  
     	Nobs(ii,jj) = tmp1[0] ;  
     	Nobs(jj,ii)=Nobs(ii,jj) ;  
     	Ntot_used(ii) += tmp1[0] ;  
     	Ntot_used(jj) += tmp1[0] ;	  
     	tmp2 = res["thresh1"] ;  
     	for (int uu=0;uu<maxK3 ;uu++){  
     		thresh(ii,uu) += tmp2[uu]*Nobs(ii,jj) ;  
     				}  
     	tmp2 = res["thresh2"] ;  
     	for (int uu=0;uu<maxK3 ;uu++){  
     		thresh(jj,uu) += tmp2[uu]*Nobs(ii,jj) ;  
     				}  
     	tmp3 = res["maxK1"];  
     	if (tmp3[0] > maxcat[ii] ){ maxcat[ii] = tmp3[0] ; }  
     	tmp3 = res["maxK2"];  
     	if (tmp3[0] > maxcat[jj] ){ maxcat[jj] = tmp3[0] ; }  
     		}  
     	}  
       
     for (int ii=0 ; ii < I;ii++){	  
     for (int uu=0;uu<maxK3 ; uu++){  
     	thresh(ii,uu) = thresh(ii,uu) / Ntot_used(ii) ;  
     	if ( uu > maxcat[ii]  ){  
     		thresh(ii,uu) = 99 ;  
     				}  
     			}  
     	rho(ii,ii) = 1 ;  
     			}  
       
     //*************************************************      
     // OUTPUT              
                   
     return Rcpp::List::create(    
         Rcpp::_["rho"] = rho  , 
         Rcpp::_["thresh"]=thresh ,  
         Rcpp::_["maxcat"] = maxcat , 
         Rcpp::_["Nobs"]=Nobs , 
         Rcpp::_["Ntot_used"] = Ntot_used  
         ) ;    
}



