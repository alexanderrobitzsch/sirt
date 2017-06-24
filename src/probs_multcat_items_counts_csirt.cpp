

#include <Rcpp.h>

using namespace Rcpp;



///********************************************************************
///** probs_pcm_groups_C
// [[Rcpp::export]]
Rcpp::List probs_pcm_groups_C( Rcpp::NumericMatrix dat, 
	Rcpp::NumericMatrix dat_resp, Rcpp::NumericVector group, 
	Rcpp::NumericMatrix probs, 
	int CC, int TP ){
 
     //*** input probs  
     // probs[ items , categories , theta points , groups ]  
       
     int I = dat.ncol();  // number of items  
     int N = dat.nrow() ; // number of subjects  
       
     // create likelihood object  
     Rcpp::NumericMatrix fyiqk(N,TP);  
     fyiqk.fill(1);    
       
     // int nn = 0 ;  
     for (int nn=0;nn<N;nn++){  
     for (int ii=0;ii<I;ii++){  
     if (dat_resp(nn,ii)==1){  
        for (int tt=0;tt<TP;tt++){  
     	// probs ( ii , cc , tt , gg ) =   
     	// probs_C(ii ,  cc + tt*CC + gg * CC*TP )  
     	fyiqk(nn,tt) = fyiqk(nn,tt) * probs( ii , dat(nn,ii) + tt * CC + group[nn] * CC * TP ) ;  
     			}  // end tt  
     		}  // end if dat_resp = 1  
     	} // end ii  
     }   // end nn  
     		  
     ///////////////////////////////////////  
     /// OUTPUT                  
       
     return Rcpp::List::create(  
     		Rcpp::_["fyiqk"] = fyiqk  
     			) ;  
}




//*******************************************
// individual likelihood (no groups)


///********************************************************************
///** probs_pcm_nogroups_C
// [[Rcpp::export]]
Rcpp::List probs_pcm_nogroups_C( Rcpp::NumericMatrix dat, 
	Rcpp::NumericMatrix dat_resp, Rcpp::NumericMatrix probs, 
	int CC, int TP ){

     //*** input probs  
     // probs[ items , categories , theta points , groups ]  
       
     int I = dat.ncol();  // number of items  
     int N = dat.nrow() ; // number of subjects  
       
     // create likelihood object  
     Rcpp::NumericMatrix fyiqk(N,TP);  
     fyiqk.fill(1);    
       
     // int nn = 0 ;  
     for (int nn=0;nn<N;nn++){  
     for (int ii=0;ii<I;ii++){  
     if (dat_resp(nn,ii)==1){  
        for (int tt=0;tt<TP;tt++){  
     	// probs ( ii , cc , tt , gg ) =   
     	// probs_C(ii ,  cc + tt*CC + gg * CC*TP )  
     	fyiqk(nn,tt) = fyiqk(nn,tt) * probs( ii , dat(nn,ii) + tt * CC  ) ;  
     			}  // end tt  
     		}  // end if dat_resp = 1  
     	} // end ii  
     }   // end nn  
     		  
     ///////////////////////////////////////  
     /// OUTPUT                                
     return Rcpp::List::create(  
     		Rcpp::_["fyiqk"] = fyiqk  
     			) ;  
}




//********************************************************
// calculation of expected counts

///********************************************************************
///** calccounts_pcm_groups_C
// [[Rcpp::export]]
Rcpp::List calccounts_pcm_groups_C( Rcpp::NumericMatrix dat, 
	Rcpp::NumericMatrix dat_resp, Rcpp::NumericVector group, 
	Rcpp::NumericMatrix fyiqk, Rcpp::NumericMatrix pik, 
	int CC, Rcpp::NumericVector weights ){
  
      
     int TP = fyiqk.ncol() ;  
     int G = pik.ncol() ;  
     int N = dat.nrow() ;  
     int I = dat.ncol() ;  
     Rcpp::NumericMatrix fqkyi(N,TP) ;  
     Rcpp::NumericMatrix pik1(TP,G);  
     double t1 = 0 ;  
       
     //**************  
     // calculate posterior  
     for (int nn=0;nn<N;nn++){  
     	t1 = 0 ;  
     	for (int tt=0;tt<TP;tt++){  
     		fqkyi(nn,tt) = fyiqk(nn,tt) * pik( tt , group[nn] ) ;  
     		t1 += fqkyi(nn,tt) ;  
     			}  
     	fqkyi.row(nn) = fqkyi.row(nn) / t1 ;  
     }  
       
     ///********  
     // calculate expected counts  
     // probs ( ii , cc , tt , gg ) =   
     // probs_C(ii ,  cc + tt*CC + gg * CC*TP )  
     Rcpp::NumericMatrix nik(I,CC*TP*G);  
       
     for (int ii=0;ii<I;ii++){  
     for (int nn=0;nn<N;nn++){  
     //int nn = 0 ;  
     if (dat_resp(nn,ii) ==1){  
     for(int tt=0;tt<TP;tt++){           
     	nik(ii,dat(nn,ii)+tt*CC+group[nn]*CC*TP ) +=  fqkyi(nn,tt) * weights[nn] ;  
     		           } // end tt           
          		}  // end if dat_resp == 1   
     		} // end nn  
         }  // end ii  
       
           
     //*****************    
     // calculate loglikelihood and updated group distributions  
       
     double LL=0;    
             
     //        ll[gg] <- sum( dat1[group==gg,2] * log( rowSums( f.yi.qk[group==gg,] *     
     //					outer( rep(1,nrow(f.yi.qk[group==gg,])) , pi.k[,gg] ) ) ) )    
              
     for (int nn=0;nn<N;++nn){    
         double total = 0; 	    
         for (int tt=0;tt<TP;++tt){    
                total += fyiqk(nn,tt) * pik(tt,group[nn]);    
                pik1(tt,group[nn]) += fqkyi(nn,tt)*weights[nn] ;   
                           } // end tt    
          LL += log( total ) * weights[nn] ;    
                      }  // end nn        
           
     ///////////////////////////////////////  
     /// OUTPUT                  
     return Rcpp::List::create(  
     		Rcpp::_["LL"] = LL ,  
     		Rcpp::_["fqkyi"] = fqkyi ,  
     		Rcpp::_["nik"] = nik ,  
                Rcpp::_["count_pik"] = pik1  
    			) ;  
}




