//// File Name: rm_smirt_mml2_code.cpp
//// File Version: 5.15
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// user includes


//**********************************************************************
////////////////////////////////////////////////////////////////////////
/// array multiplication
/// The following function is used in arraymult1 in file rm.hrm_alg
////////////////////////////////////////////////////////////////////////
//**********************************************************************



///********************************************************************
///** probs_pcm_nogroups_C
// [[Rcpp::export]]
Rcpp::NumericMatrix rm_arraymult1( Rcpp::NumericMatrix AA, 
	Rcpp::NumericVector dimAA, Rcpp::NumericMatrix BB, 
	Rcpp::NumericVector dimBB ){
           
      Rcpp::NumericMatrix CC( dimAA[0]*dimAA[1] , dimBB[2] ) ; 
      int a1 = dimAA[0];
      int a2 = dimAA[1];
      int a3 = dimAA[2];
      int b3 = dimBB[2]; 
      
      // ii -> loop within a matrix
      //for (int ii=0 ; ii < a1 ; ++ii ){
      for (int zz=0 ; zz < a2 ; ++zz){
          for (int ii=0 ; ii < a1 ; ++ii ){        
              for (int hh=0 ; hh < b3 ; ++hh){ // loop over columns
                  for (int kk=0 ; kk < a3 ; ++kk){
                      CC(ii+zz*a1,hh) += AA(ii+zz*a1,kk)*BB(ii+a1*kk,hh)  ; // *BB(kk,hh) ;
                              }
                          }
                      }
                  }
      
      // output
      return CC;
}



//**********************************************************************
////////////////////////////////////////////////////////////////////////
// calculation posterior distribution
// RM_CALCPOST used in R function rm_calclike in file rm.alg
////////////////////////////////////////////////////////////////////////
//**********************************************************************



///********************************************************************
///** RM_CALCPOST
// [[Rcpp::export]]
Rcpp::List RM_CALCPOST( Rcpp::NumericMatrix DAT2, 
	Rcpp::NumericMatrix DAT2RESP, Rcpp::NumericMatrix PROBS, 
	Rcpp::NumericVector KK ){
         
     int N=DAT2.nrow();  
     int I=DAT2.ncol();  
     int TP=PROBS.ncol();  
     int KKK = KK[0] + 1 ;   
       
     //*****  
     // calculate individual likelihood  
     Rcpp::NumericMatrix fyiqk (N,TP) ;  
     fyiqk.fill(1);  
     for (int ii=0;ii<I;++ii){      
     for (int nn=0;nn<N;++nn){  
         if ( DAT2RESP(nn,ii)>0){  
         for (int tt=0;tt<TP;++tt){  
             fyiqk(nn,tt) = fyiqk(nn,tt) * PROBS( KKK*ii + DAT2(nn,ii) , tt ) ;  
                         }  
                     }  
                 }  
             }  
       		  
     ///////////////////////////////////////////////////////  
     ///////////// O U T P U T   ///////////////////////////  
     return Rcpp::List::create(
            Rcpp::_["fyiqk"] = fyiqk 
                );  		
}





//**********************************************************************
////////////////////////////////////////////////////////////////////////
/// calculation of probabilities
// probraterfct1 in file rm.hrm_alg
////////////////////////////////////////////////////////////////////////
//**********************************************************************


///********************************************************************
///** rm_probraterfct1
// [[Rcpp::export]]
Rcpp::List rm_probraterfct1( Rcpp::NumericMatrix CRA, 
	Rcpp::NumericVector DRA , Rcpp::IntegerVector dimAA, 
	Rcpp::NumericMatrix BB, Rcpp::IntegerVector dimBB ){

      int K = CRA.ncol();
      int I = CRA.nrow();
      
      // define input and arrays
      // dimA (I,K+1,K+1)
      // Rcpp::IntegerVector dimAA(dimA);
      // array B (I,K+1,TP)
      // Rcpp::NumericMatrix BB(B);
      // Rcpp::IntegerVector dimBB(dimB);
      Rcpp::NumericMatrix CC( dimAA[0]*dimAA[1] , dimBB[2] ) ;
      
      ////////////////////////////////////////////
      
      //*********************
      // create h1 matrix
      Rcpp::NumericMatrix h1(I*(K+1),K) ;
      int nrh1 = h1.nrow();
      
      for (int kk=0;kk<(K+1);++kk){
          for (int ii=0;ii<I;++ii){
              for (int jj=0;jj<K;++jj){
                  h1(ii+I*kk,jj) = CRA( ii , jj )  - DRA[ii]*kk ;
                                     }
                                 }
                          }
                          
      //**********************************
      // compute logistic distribution
      // for (int kk=0; kk<K; ++kk){
      for (int kk=0;kk<K;++kk){
          Rcpp::NumericVector inpv =h1(_,kk) ;
          Rcpp::NumericVector res0=plogis(inpv) ;
          for (int ii=0;ii<nrh1;++ii){
              h1(ii,kk) = res0[ii] ;
                               }
                      }
      
      //***********************************
      // compute matrix with rater probabilities  
      
      // P(X=0) ;                
      Rcpp::NumericMatrix PRA(I*(K+1),K+1) ;         
      for (int jj=0;jj<(K+1);++jj){                
          for (int ii=0;ii<I;++ii){
              PRA(ii,jj) = h1(ii+I*jj, 0 ) ;
                      }    
                  }
      // other categories
      for (int cc=1;cc<K;++cc){            
      for (int jj=0;jj<(K+1);++jj){            
          for (int ii=0;ii<I;++ii){
              PRA(ii+I*cc,jj) = h1(ii+I*jj, cc ) - h1(ii+I*jj, cc-1 ) ;
                      }               
                  }
              }
      // last category
      int cc=K ;        
      for (int jj=0;jj<(K+1);++jj){            
          for (int ii=0;ii<I;++ii){
              PRA(ii+I*cc,jj) = 1 - h1(ii+I*jj, cc-1 ) ;
                      }               
                  }
      
      //**************************
      // multiplication 
      Rcpp::NumericMatrix AA=PRA ;                   
      int a1 = dimAA[0];
      int a2 = dimAA[1];
      int a3 = dimAA[2];
      int b3 = dimBB[2]; 
      
      // ii -> loop within a matrix
      //for (int ii=0 ; ii < a1 ; ++ii ){
      for (int zz=0 ; zz < a2 ; ++zz){
          for (int ii=0 ; ii < a1 ; ++ii ){        
              for (int hh=0 ; hh < b3 ; ++hh){ // loop over columns
                  for (int kk=0 ; kk < a3 ; ++kk){
                      CC(ii+zz*a1,hh) += AA(ii+zz*a1,kk)*BB(ii+a1*kk,hh)  ; // *BB(kk,hh) ;
                              }
                          }
                      }
                  }
                           
      ///////////////////////////////////////////////////////
      ///////////// O U T P U T   ///////////////////////////
      return Rcpp::List::create(
         Rcpp::_["PRA"] = PRA , 
         Rcpp::_["probtotal"] = CC
                    );            
                  
}






//********************************************************************
//********************************************************************
//********************************************************************
// rm_facets_calcprobs
//********************************************************************
//********************************************************************
//********************************************************************

///********************************************************************
///** rm_facets_calcprobs_cpp
// [[Rcpp::export]]
Rcpp::NumericMatrix rm_facets_calcprobs_cpp( Rcpp::NumericVector b_item, 
	Rcpp::NumericVector b_rater, Rcpp::NumericMatrix Qmatrix, 
	Rcpp::NumericMatrix tau_item, int K, int I, int TP, 
	Rcpp::NumericVector a_item, Rcpp::NumericVector a_rater, 
	Rcpp::NumericVector item_index, Rcpp::NumericVector rater_index, 
	Rcpp::NumericVector theta_k ){

     // int RR = as<int>(RR_);  
     //    probs <- .rm.facets.calcprobs( b.item , b.rater , Qmatrix , tau.item ,  
     //           VV , K , I , TP , a.item , a.rater , item.index , rater.index ,  
     //           theta.k ,RR )       
       
     //***** calculate b  
     // b <- tau.item[ item.index , ]  
     Rcpp::NumericMatrix b(I,K) ;  
     // b0 <- ( matrix( b.rater , nrow= RR , ncol=K) )[ rater.index , ] * 	Qmatrix[ item.index ,]	   
     // b <- b + b0  
     for (int ii=0; ii<I ; ii++){  
        b.row(ii) = tau_item.row( item_index[ii] ) ;  
        for (int kk=0;kk<K;kk++){   	     
         b(ii,kk) = b(ii,kk) + b_rater[ rater_index[ii] ] * Qmatrix( item_index[ii] , kk ) ;  
         			}  
         		}  
       
     //****** calculate a    		  
     // a <- a.item[ item.index ] * a.rater[ rater.index ]    		  
     Rcpp::NumericVector a(I) ;  
     for (int ii=0;ii<I;ii++){  
     	a[ii] = a_item[ item_index[ii] ] * a_rater[ rater_index[ii] ] ;  
     			}  
     //******* calculate modified Q-matrix  
     // Qmatrix=Qmatrix[item.index,]  
     Rcpp::NumericMatrix Q(I,K) ;  
     for (int ii=0;ii<I;ii++){  
     	Q.row(ii) = Qmatrix.row( item_index[ii] ) ;  
     			}  
       
     ///**************************  
     // compute response probabilities according to the generalized partial credit model  
       
     //     probs <- array( 0 , dim=c(I,K+1,TP) )   # categories 0 , ... , K  
     Rcpp::NumericMatrix probs(I,(K+1)*TP) ;  
     //    for (kk in 1:K){  
     //        l0 <- matrix( - b[,kk] , nrow=I,ncol=TP)  
     //        l0 <- l0 + outer( a * Qmatrix[ , kk] , theta.k )  
     //        probs[,kk+1,] <- l0  
     //                }  
     // probs(ii , kk , tt ) ~ probs( ii , kk + tt * (K+1) )  
     double tmp1 = 0 ;   
       
     for (int tt=0;tt<TP;tt++){  
     for (int ii=0;ii<I;ii++){  
     	tmp1 = 1 ;	  
     	probs( ii , tt*(K+1) ) = 1 ;  
     	for (int kk=0;kk<K;kk++){  
     		probs( ii , kk+1 + tt*(K+1) ) = exp( - b(ii,kk) + a[ii] * Q(ii,kk) * theta_k[tt] );  
     		tmp1 += probs( ii , kk+1 + tt*(K+1) ) ;  
     				}  
     	for (int kk=0;kk<K+1;kk++){  
     		probs( ii , kk + tt*(K+1) ) = probs( ii , kk + tt*(K+1) ) / tmp1 ;  
     				}								  
     		}   // end ii  
     	} // end tt  
       
     return probs ;       	  
}






//**********************************************************************
////////////////////////////////////////////////////////////////////////
/// SMIRT_CALCPOST
/// smirt - calculation of posterior probabilities
////////////////////////////////////////////////////////////////////////
//**********************************************************************

///********************************************************************
///** SMIRT_CALCPOST
// [[Rcpp::export]]
Rcpp::List SMIRT_CALCPOST( Rcpp::IntegerMatrix DAT2, 
	Rcpp::IntegerMatrix DAT2RESP, Rcpp::NumericMatrix PROBS, 
	Rcpp::NumericMatrix DAT2IND, Rcpp::NumericVector PIK, 
	Rcpp::NumericVector KK1 ){
         
     int N=DAT2.nrow();  
     int I=DAT2.ncol();  
     int TP=PROBS.ncol();  
     int KK=KK1[0] ;   
       
     //*****  
     // calculate individual likelihood  
     Rcpp::NumericMatrix fyiqk (N,TP) ;  
     fyiqk.fill(1);  
     for (int ii=0;ii<I;++ii){      
     for (int nn=0;nn<N;++nn){  
         if ( DAT2RESP(nn,ii)>0){  
         for (int tt=0;tt<TP;++tt){  
             fyiqk(nn,tt) = fyiqk(nn,tt) * PROBS( 2*ii + DAT2(nn,ii) , tt ) ;  
                         }  
                     }  
                 }  
             }  
       
     //****  
     // calculate posterior	  
     Rcpp::NumericMatrix fqkyi (N,TP) ;  
     for (int nn=0;nn<N;++nn){  
     	double total = 0; 	  
     	for (int tt=0;tt<TP;++tt){  
     		fqkyi(nn,tt) = fyiqk(nn,tt)*PIK[tt] ;		  
     		total += fqkyi(nn,tt) ;  
     			}  
     	for (int tt=0;tt<TP;++tt){  
     		fqkyi(nn,tt) = fqkyi(nn,tt)/total ;  
     			}  
     		}  
     //*****  
     // calculate counts  
     for (int tt=0;tt<TP ;++tt){  
     	PIK[tt] = 0 ;  
     	for (int nn=0;nn<N;++nn){	  
     		PIK[tt] += fqkyi(nn,tt) ;  
     				}  
     	PIK[tt] = PIK[tt] / N ;   
     			}  
       
     Rcpp::NumericMatrix nik (TP, I*(KK+1)) ;  
     Rcpp::NumericMatrix NIK (TP, I) ;			  
     for (int tt=0;tt<TP;++tt){  
     for (int ii=0; ii < I ; ++ii){  
     for (int kk=0;kk<KK+1;++kk){			  
     	for (int nn=0;nn<N;++nn){  
     //			nik( tt , ii + kk*I  ) += DAT2RESP(nn,ii)*(DAT2(nn,ii)==kk)*fqkyi(nn,tt)  ;  
     			nik( tt , ii + kk*I  ) += DAT2IND(nn,ii+kk*I) *fqkyi(nn,tt)  ;  
     				}  // end nn  
     		NIK(tt,ii) += nik(tt,ii+kk*I ) ; 				  
     			} // end kk  
     		}  // end ii  
     	}  // end tt  
     			  
     ///////////////////////////////////////////////////////  
     ///////////// O U T P U T   ///////////////////////////  
     return Rcpp::List::create(
        Rcpp::_["fyiqk"] = fyiqk , 
        Rcpp::_["f.qk.yi"]=fqkyi ,   
     	Rcpp::_["pi.k"] = PIK , 
        Rcpp::_["n.ik"] = nik , 
        Rcpp::_["N.ik"]=NIK 
                 );   
     	  
}



//**********************************************************************
////////////////////////////////////////////////////////////////////////
/// SMIRT_CALCPROB_COMP
/// smirt - calculation of probabilities in the compensatory model
////////////////////////////////////////////////////////////////////////
//**********************************************************************

///********************************************************************
///** SMIRT_CALCPROB_COMP
// [[Rcpp::export]]
Rcpp::NumericMatrix SMIRT_CALCPROB_COMP( Rcpp::NumericMatrix A, 
	Rcpp::NumericVector B, Rcpp::NumericMatrix QQ, 
	Rcpp::NumericMatrix THETA, Rcpp::NumericVector CC, 
	Rcpp::NumericVector DD ){
 
    int I=A.nrow();
    int D=A.ncol();
    int TP=THETA.nrow();

    // create matrix of probabilities
    Rcpp::NumericMatrix prob (I,TP) ;
    prob.fill(1);
    Rcpp::NumericVector yy (TP);

    for (int ii=0;ii<I;++ii){
            for (int tt=0; tt<TP; ++tt){
            yy[tt] = - B[ii] ; 	
           for (int dd=0;dd<D;++dd){
                   yy[tt] += A(ii,dd)*QQ(ii,dd)*THETA(tt,dd)  ;
                    } // end dd 
                        }   // end tt
            Rcpp::NumericVector p1=plogis(yy) ;
            for (int tt=0; tt<TP; ++tt){            
                prob(ii,tt)= p1[tt] ;
                  }  // end tt
    //***
    // include guessing and slipping parameters
     if ( ( CC[ii] > 0 ) || ( DD[ii] < 1 ) ){ 
        for (int tt=0;tt<TP;++tt){    
            prob(ii,tt) = CC[ii] + ( DD[ii]-CC[ii] )* prob(ii,tt) ;
                    }   // end tt
                } // end if condition for guessing or slipping                
        }       // end ii
        
    ///////////////////////////////////////
    /// OUTPUT    
    return prob ;             
}



//**********************************************************************
////////////////////////////////////////////////////////////////////////
/// SMIRT_CALCPROB_NONCOMP
/// smirt - calculations in the noncompensatory model
////////////////////////////////////////////////////////////////////////
//**********************************************************************


///********************************************************************
///** SMIRT_CALCPROB_NONCOMP
// [[Rcpp::export]]
Rcpp::NumericMatrix SMIRT_CALCPROB_NONCOMP( Rcpp::NumericMatrix A, 
	Rcpp::NumericMatrix B, Rcpp::NumericMatrix QQ, 
	Rcpp::NumericMatrix THETA, Rcpp::NumericVector CC, 
	Rcpp::NumericVector DD ){
        
     int I=A.nrow();  
     int D=A.ncol();  
     int TP=THETA.nrow();  
       
     // create matrix of probabilities  
     Rcpp::NumericMatrix prob (I,TP) ;  
     prob.fill(1);  
     Rcpp::NumericVector yy (TP);  
       
     for (int ii=0;ii<I;++ii){  
     for (int dd=0;dd<D;++dd){  
         if ( QQ(ii,dd)>0 ){  
             for (int tt=0; tt<TP; ++tt){  
                    yy[tt] =  A(ii,dd)*QQ(ii,dd)*THETA(tt,dd)-B(ii,dd)   ;  
                         }   // end tt  
             Rcpp::NumericVector p1=plogis(yy) ;  
             for (int tt=0; tt<TP; ++tt){              
                 prob(ii,tt)=prob(ii,tt)*p1[tt] ;  
                   }  // end tt  
                 }       // end if Q(ii,dd)>0  
             }          // end dd   
     //***  
     // include guessing and slipping parameters  
      if ( ( CC[ii] > 0 ) || ( DD[ii] < 1 ) ){   
         for (int tt=0;tt<TP;++tt){      
             prob(ii,tt) = CC[ii] + ( DD[ii]-CC[ii] )* prob(ii,tt) ;  
                     }   // end tt  
                 } // end if condition for guessing or slipping                  
         }       // end ii  
           
     ///////////////////////////////////////  
     /// OUTPUT      
     return prob ;               
}



//**********************************************************************
////////////////////////////////////////////////////////////////////////
/// SMIRT_CALCPROB_PARTCOMP
/// smirt - partially noncompensatory model
////////////////////////////////////////////////////////////////////////
//**********************************************************************

///********************************************************************
///** SMIRT_CALCPROB_PARTCOMP
// [[Rcpp::export]]
Rcpp::NumericMatrix SMIRT_CALCPROB_PARTCOMP( Rcpp::NumericMatrix A, 
	Rcpp::NumericMatrix B, Rcpp::NumericMatrix QQ, 
	Rcpp::NumericMatrix THETA, Rcpp::NumericVector CC, 
	Rcpp::NumericVector DD , Rcpp::NumericVector MUI ){
       
     int I=A.nrow();  
     int D=A.ncol();  
     int TP=THETA.nrow();         
     // create matrix of probabilities  
     Rcpp::NumericMatrix prob (I,TP) ;  
     prob.fill(1);  
     Rcpp::NumericVector yy1 (1);  
     Rcpp::NumericVector yy2 (1); 
     Rcpp::NumericVector yy3 (1);
     double tmp1 ;
     
     for (int ii=0;ii<I;++ii){  
     for (int tt=0; tt<TP; ++tt){
     	    yy1[0] = 0 ; 
            yy2[0] = 1 ;    	    
     for (int dd=0;dd<D;++dd){  
         if ( QQ(ii,dd)>0 ){    
             tmp1 = A(ii,dd)*QQ(ii,dd)*THETA(tt,dd)-B(ii,dd) ;
             yy1[0] = yy1[0] + tmp1  ;
             yy2[0] = yy2[0] * ( 1 + exp( tmp1 ) ) ;
                 }       // end if Q(ii,dd)>0  
             yy3[0] = exp( yy1[0] ) ;                              
                 }          // end dd   
          prob(ii,tt) = yy3[0]/(  MUI[ii] * yy2[0] + (1-MUI[ii])*(1+yy3[0]) ) ;
             
     //***  
     // include guessing and slipping parameters  
      if ( ( CC[ii] > 0 ) || ( DD[ii] < 1 ) ){        
             prob(ii,tt) = CC[ii] + ( DD[ii]-CC[ii] )* prob(ii,tt) ;  
                 } // end if condition for guessing or slipping                   
     		} // end tt                  
         }       // end ii  
           
     ///////////////////////////////////////  
     /// OUTPUT      
     return prob ;               
}



/////////////////////////////////////////////////////////////
/// rasch.mml2 : calculate counts

///********************************************************************
///** MML2_RASCHTYPE_COUNTS
// [[Rcpp::export]]
Rcpp::List MML2_RASCHTYPE_COUNTS( Rcpp::NumericMatrix DAT2, 
	Rcpp::NumericMatrix DAT2RESP, Rcpp::NumericVector DAT1, 
	Rcpp::NumericMatrix FQKYI, Rcpp::NumericVector PIK, 
	Rcpp::NumericMatrix FYIQK ){

     int N=DAT2.nrow();  
     int I=DAT2.ncol();  
     int TP=FQKYI.ncol();  
       
     //**********************************  
     // calculate total counts n.k  
     Rcpp::NumericVector NK (TP);  
     NK.fill(0);  
     for (int tt=0;tt<TP;++tt){      
         for (int nn=0;nn<N;++nn){  
     //        NK[tt] = NK[tt] + FQKYI(nn,tt)*DAT1[nn] ;   
               NK[tt] += FQKYI(nn,tt)*DAT1[nn] ;   
                         } // end nn  
                     }  // end tt  
       
     //**********************************  
     // calculate n.jk and r.jk  
     Rcpp::NumericMatrix NJK (I,TP) ;  
     Rcpp::NumericMatrix RJK (I,TP) ;  
     NJK.fill(0) ;                   
     RJK.fill(0) ;    
                       
     for (int ii=0;ii<I;++ii){  
     for (int tt=0;tt<TP;++tt){  
         for (int nn=0;nn<N;++nn){  
             if (DAT2RESP(nn,ii)>0){  
                 NJK(ii,tt) += DAT1[nn] * FQKYI(nn,tt) ;   
                 RJK(ii,tt) += DAT1[nn] * FQKYI(nn,tt) * DAT2(nn,ii) ;               
                             }       // end if ...  
             }           //     end nn  
         } // end tt  
     }    // end ii  
       
     //*****************  
     // calculate loglikelihood  
     double LL=0;  
     //        ll[gg] <- sum( dat1[group==gg,2] * log( rowSums( f.yi.qk[group==gg,] *   
     //					outer( rep(1,nrow(f.yi.qk[group==gg,])) , pi.k[,gg] ) ) ) )  
     for (int nn=0;nn<N;++nn){  
         double total = 0; 	  
         for (int tt=0;tt<TP;++tt){  
               total += FYIQK(nn,tt) * PIK[tt] ;  
                             } // end tt  
         LL += log( total )*DAT1[nn] ;  
                 }  // end nn  
                                                 
     ///////////////////////////////////////////////////////  
     ///////////// O U T P U T   ///////////////////////////  
     return Rcpp::List::create(
         Rcpp::_["nk"] = NK , 
         Rcpp::_["njk"] = NJK ,   
         Rcpp::_["rjk"] = RJK , 
         Rcpp::_["ll"] = LL 
                 );   
}




//////////////////////////////////////////////////////////////////
/// rasch.mml2 - calculation of posterior distribution

///********************************************************************
///** MML2_CALCPOST_V1
// [[Rcpp::export]]
Rcpp::List MML2_CALCPOST_V1( Rcpp::NumericMatrix DAT2, 
	Rcpp::NumericMatrix DAT2RESP, Rcpp::NumericMatrix PROBS ){
       
     int N=DAT2.nrow();  
     int I=DAT2.ncol();  
     int TP=PROBS.ncol();  
       
     //*****  
     // calculate individual likelihood  
     Rcpp::NumericMatrix fyiqk (N,TP) ;  
     fyiqk.fill(1);  
     for (int ii=0;ii<I;++ii){      
     for (int nn=0;nn<N;++nn){  
         if ( DAT2RESP(nn,ii)>0){  
         for (int tt=0;tt<TP;++tt){  
             fyiqk(nn,tt) = fyiqk(nn,tt) * PROBS( 2*ii + DAT2(nn,ii) , tt ) ;                         
                         }  
                     }  
                 }  
             }  
        			  
     ///////////////////////////////////////////////////////  
     ///////////// O U T P U T   ///////////////////////////  
     return Rcpp::List::create(
       Rcpp::_["fyiqk"] = fyiqk 
                   );   

}




//////////////////////////////////////////////////////////////////
/// rasch.mml2 - calculation of posterior distribution (V2)

///********************************************************************
///** MML2_CALCPOST_V2
// [[Rcpp::export]]
Rcpp::List MML2_CALCPOST_V2( Rcpp::NumericMatrix DAT2, 
	Rcpp::NumericMatrix DAT2RESP, Rcpp::NumericMatrix PROBS ){
       
     int N=DAT2.nrow();  
     int I=DAT2.ncol();  
     int TP=PROBS.ncol();  
     //*****  
     // calculate individual likelihood  
     Rcpp::NumericMatrix fyiqk (N,TP) ;  
     fyiqk.fill(1);  
     for (int ii=0;ii<I;++ii){      
     for (int nn=0;nn<N;++nn){  
         if ( DAT2RESP(nn,ii)>0){  
         for (int tt=0;tt<TP;++tt){  
//             fyiqk(nn,tt) = fyiqk(nn,tt) * PROBS( 2*ii + DAT2(nn,ii) , tt ) ;
	if ( ( DAT2(nn,ii) < 1 ) & ( DAT2(nn,ii) > 0 ) ){
             fyiqk(nn,tt) = fyiqk(nn,tt) * pow( PROBS(2*ii+1,tt),DAT2(nn,ii) )*
                	pow( PROBS( 2*ii  , tt ) , 1 - DAT2(nn,ii) ) ;
         } else {
             fyiqk(nn,tt) = fyiqk(nn,tt) * PROBS( 2*ii + DAT2(nn,ii) , tt ) ;                			
                		}
                         }  
                     }  
                 }  
             }  
         
     ///////////////////////////////////////////////////////  
     ///////////// O U T P U T   ///////////////////////////  
     return Rcpp::List::create(
        Rcpp::_["fyiqk"] = fyiqk 
           );   				
}




//////////////////////////////////////////////////////////////////
/// rasch.mml2 - calculation of posterior distribution 
/// fuzzy log likelihood in the belief function framework


///********************************************************************
///** MML2_CALCPOST_V3
// [[Rcpp::export]]
Rcpp::List MML2_CALCPOST_V3( Rcpp::NumericMatrix DAT2,
	Rcpp::NumericMatrix DAT2RESP, Rcpp::NumericMatrix PROBS ){
       
     int N=DAT2.nrow();  
     int I=DAT2.ncol();  
     int TP=PROBS.ncol();       
     //*****  
     // calculate individual likelihood  
     Rcpp::NumericMatrix fyiqk (N,TP) ;  
     fyiqk.fill(1);  
     for (int ii=0;ii<I;++ii){      
     for (int nn=0;nn<N;++nn){  
         if ( DAT2RESP(nn,ii)>0){  
         for (int tt=0;tt<TP;++tt){  
//             fyiqk(nn,tt) = fyiqk(nn,tt) * PROBS( 2*ii + DAT2(nn,ii) , tt ) ;
	if ( ( DAT2(nn,ii) < 1 ) & ( DAT2(nn,ii) > 0 ) ){
             fyiqk(nn,tt) = fyiqk(nn,tt) * 
              (PROBS(2*ii+1,tt) * DAT2(nn,ii) + (PROBS( 2*ii , tt ) *(1-DAT2(nn,ii))) );
                		} else {
             fyiqk(nn,tt) = fyiqk(nn,tt) * PROBS( 2*ii + DAT2(nn,ii) , tt ) ;                			
                		}
                         }  
                     }  
                 }  
             }  
            			  
     ///////////////////////////////////////////////////////  
     ///////////// O U T P U T   ///////////////////////////  
     return Rcpp::List::create(
         Rcpp::_["fyiqk"] = fyiqk 
                        );   
     	  
}









