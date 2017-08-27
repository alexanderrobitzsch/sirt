//// File Name: gooijer_isop.cpp
//// File Version: 5.07
//// File Last Change: 2017-02-18 18:49:46


// [[Rcpp::depends(RcppArmadillo)]]


// includes from the plugin
#include <RcppArmadillo.h>
#include <Rcpp.h>


// [[packageincludes]]
#include "gooijercsntableaux.h"
// #include "c:/Users/robitzsch/Dropbox/Eigene_Projekte/R-Routinen/IRT-Functions/sirt_package/1.15/sirt_work/src/gooijercsntableaux__1.08.h"



using namespace Rcpp;


//********************************************************************
//********************************************************************
//********************************************************************
// unidimensionality test of Gooijer
//********************************************************************
//********************************************************************
//********************************************************************


///********************************************************************
///** gooijer_csn_table
// [[Rcpp::export]]
Rcpp::List gooijer_csn_table( Rcpp::NumericMatrix dat, 
	Rcpp::NumericMatrix dat_perm, int RR , int NS, 
	int progress, Rcpp::NumericVector progress_vec, 
	Rcpp::NumericMatrix score_index ){
 
     int N=dat.nrow() ;  
     Rcpp::NumericVector stat_perm(RR) ;  
     Rcpp::NumericMatrix sampleM = dat_perm ;  
     // compute Gooijer statistic for original data  
     Rcpp::NumericVector stat = gta(dat) ;  
     //******* ;  
     Rcpp::NumericVector s1(1);  
     int zz=0;  
     if ( progress==1){  
     	Rcpp::Rcout << "|" <<   std::flush ;  
     			}  
     for (int rr=0;rr<RR;rr++){  
     	// permutation of original data ;  
     	sampleM = gooijer_permutation(sampleM , NS , N , score_index ) ;  
     	// compute statistic for permuted data ;  
     	s1 = gta(sampleM ) ;  
     	stat_perm[rr] = s1[0] ;  
     	if ( (progress==1) & ( rr==progress_vec[zz] ) ){  
     		zz = zz+1 ;  
     		if ( zz==10){	zz = 9 ; }   
     		Rcpp::Rcout << "-" <<   std::flush ;  
     				}  
     		}  
     if ( progress==1){  
     	Rcpp::Rcout << "|" <<   std::flush << std::endl ;  
     			}			  
      		  
        		  
     //*************************************************      
     // OUTPUT                                 
      return Rcpp::List::create(    
         Rcpp::_["stat"] = stat ,    
         Rcpp::_["stat_perm"]=stat_perm  
         ) ;  
}
//********************************************************************
//********************************************************************
//********************************************************************




//********************************************************************
//********************************************************************
//********************************************************************
// ISOP test
//********************************************************************
//********************************************************************
//********************************************************************



///********************************************************************
///** isop_tests_C
// [[Rcpp::export]]
Rcpp::List isop_tests_C( Rcpp::NumericMatrix dat, 
	Rcpp::NumericMatrix dat_resp, Rcpp::NumericVector weights, 
	Rcpp::NumericVector jackunits, int JJ ){
  
       
     int N=dat.nrow();  
     int I=dat.ncol();  
       
     Rcpp::NumericMatrix Esi(I,JJ+1);  
     Rcpp::NumericMatrix Edi(I,JJ+1);  
     Rcpp::NumericMatrix W1i(I,JJ+1);  
     Rcpp::NumericVector W1test(JJ+1);  
       
       
     // int ii=0 ;  
     for (int ii=0;ii<I;ii++){  
     for (int nn=0;nn<N;nn++){  
     for (int mm=0;mm<N;mm++){  
     if (mm!=nn){  
     for (int jj=0;jj<I;jj++){  
         if ( jj != ii ){  
             if ( ( dat(nn,ii) > dat(mm,ii) ) & ( dat(nn,jj) > dat(mm,jj) ) ){  
     //            Esi( ii , 0 ) += weights[nn] ;  
                 Esi.row(ii) = Esi.row(ii) + weights[nn] ;  
                 Esi(ii, jackunits[nn]+1 ) +=  - weights[nn] ;  
                 if ( jackunits[nn] != jackunits[mm] ){   
                        Esi( ii , jackunits[mm] + 1 ) += - weights[nn] ;                                   
                                            }  
                                     }  
             if ( ( dat(nn,ii) < dat(mm,ii) ) & ( dat(nn,jj) < dat(mm,jj) ) ){  
                 Esi.row(ii) = Esi.row(ii) + weights[nn] ;  
                 Esi(ii, jackunits[nn]+1 ) +=  - weights[nn] ;  
                 if ( jackunits[nn] != jackunits[mm] ){   
                        Esi( ii , jackunits[mm] + 1 ) += - weights[nn] ;                                   
                                            }  
                                     }  
             if ( ( dat(nn,ii) < dat(mm,ii) ) & ( dat(nn,jj) > dat(mm,jj) ) ){  
     //            Edi( ii , 0 ) += weights[nn] ;  
                 Edi.row(ii) = Edi.row(ii) + weights[nn] ;  
                 Edi(ii, jackunits[nn]+1 ) +=  - weights[nn] ;  
                 if ( jackunits[nn] != jackunits[mm] ){   
                        Edi( ii , jackunits[mm] + 1 ) += - weights[nn] ;                                   
                                            }  
                                     }  
             if ( ( dat(nn,ii) > dat(mm,ii) ) & ( dat(nn,jj) < dat(mm,jj) ) ){  
                 Edi.row(ii) = Edi.row(ii) + weights[nn] ;  
                 Edi(ii, jackunits[nn]+1 ) +=  - weights[nn] ;  
                 if ( jackunits[nn] != jackunits[mm] ){   
                        Edi( ii , jackunits[mm] + 1 ) += - weights[nn] ;                                   
                                            }  
                                     }                                  
                                 }  // end if jj != ii  
                          }  // end for jj  
                      }  // end if mm!= nn  
                  }   // end for mm  
              } // end for nn  
          } // end for ii  
       
     // compute W1i  
     for (int ii=0;ii<I;ii++){       
         //int ii = 0 ;       
         W1i.row(ii) = ( Esi.row(ii) - Edi.row(ii)  ) /  ( Esi.row(ii) + Edi.row(ii)  ) ;  
             }  
       
     // compute statistic for the whole test          
     double tmp1 ;  
     double tmp2 ;        
     for (int jj=0;jj<JJ+1;jj++){          
         tmp2 = 0 ;  
         for (int ii=0;ii<I;ii++){          
             tmp1 =  ( Esi(ii,jj) + Edi(ii,jj) )   ;          
             tmp2 += tmp1 ;  
             W1test[jj] += tmp1 * W1i(ii,jj) ;       
                 }  
         W1test[jj] = W1test[jj] / tmp2 ;          
             }  
                              
     ///////////////////////////////////////  
     /// OUTPUT                  
     return List::create(  
             Rcpp::_["W1test"] = W1test ,   
             Rcpp::_["W1i"] = W1i ,  
             Rcpp::_["Esi"] = Esi ,  
             Rcpp::_["Edi"] = Edi  
     			) ;  
}
//********************************************************************
//********************************************************************
//********************************************************************


