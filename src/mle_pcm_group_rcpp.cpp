//// File Name: mle_pcm_group_rcpp.cpp
//// File Version: 2.06


// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


// user includes

///********************************************************************
///** mle_pcm_group_C
// [[Rcpp::export]]
Rcpp::List mle_pcm_group_C( Rcpp::NumericMatrix dat, Rcpp::NumericMatrix dat_resp, Rcpp::NumericMatrix groupM, 
	Rcpp::NumericMatrix b, Rcpp::NumericVector a, Rcpp::NumericVector maxK, Rcpp::NumericVector theta0, 
	double conv, int maxiter )
{

	int I = dat.ncol() ; 
	double eps=1e-10;
	double l1=0;  
	double l2=0;  
	double lam1=0;  
	double lam2=0;  
	double delta=0;  
	double sigdelta=0; 
	int K = b.ncol() +1;  
	int G = groupM.nrow() ;  

	Rcpp::NumericVector theta(G) ;  
	Rcpp::NumericVector setheta(G) ;  
	Rcpp::NumericVector u(I) ;  
	Rcpp::NumericVector u_resp(I) ;  
	Rcpp::NumericVector pij0(K) ;  
	Rcpp::NumericVector pij(K);  
	Rcpp::NumericVector pij2(K);  
	Rcpp::NumericVector pij3(K) ;  
	Rcpp::NumericVector Niter(G) ;  

	double th1 = 0;  
	double t1=1;  
	double absdelta = 1000 ;  
	int iter = 0 ;  

	for (int gg=0;gg<G;gg++){
		iter=0;  
		absdelta=100;  
		// iterations  
		th1 = theta0[gg] ;  
		iter=0;  
		//++++++++++++ BEGIN ITERATIONS person nn ++++++++++++++++++++++++    
		while( ( iter < maxiter ) & ( absdelta > conv)  ){  
			l1=0;
			l2=0;
			for (int nn= groupM(gg,0)-1 ; nn< groupM(gg,1) ; nn++){  
				u = dat.row(nn) ;  
				u_resp = dat_resp.row(nn) ;  
				//************** begin ITEM ii **************  
				for (int ii=0; ii<I ; ii++){  
					// calculate probabilities          
					pij3[0] = 1;  
					pij[0] = 0 ;  
					t1=1;
					for (int cc=1; cc < maxK[ii] ; cc++){  
						pij[cc] = pij[cc-1] + ( a[ii]* th1 - b(ii,cc-1) ) ;  
						pij3[cc] = exp( pij[cc] ) ;  
						t1 = pij3[cc] + t1 ;   
					}  
					for (int cc= 0 ; cc < maxK[ii] ; cc++){  
						pij[cc] = pij3[cc] / t1 ;  
					}  
					// expected value  
					lam1=0;  
					lam2=0;  
					for (int cc= 1;cc < maxK[ii];cc++){  
						lam1 = lam1 + cc * pij[cc] ;  
						lam2 = lam2 + cc*cc * pij[cc] ;  
					}  
					l1 = l1 + u_resp[ii] * a[ii] * ( u[ii] - lam1 );  
					l2 = l2 - u_resp[ii] * a[ii]*a[ii] * ( lam2 - lam1 * lam1 ) ;  
				}  // end ii  
				//************** end ITEM ii **************  
			} // end nn              
			// increment                  
			delta = - l1 / ( l2 + eps ) ;  
			sigdelta = 1 ;  // signum delta  
			if ( delta < 0 ){ sigdelta = -1 ; }  
			while( delta * sigdelta > 2){   delta = delta / 2 ; }   
				th1 = th1 + delta ;  
				absdelta = delta * sigdelta ;  
				iter ++ ; // increment iter  
		}  // end iter  
		theta[gg] = th1 ;            
		setheta[gg] = sqrt( - 1 / l2 ) ;  
		Niter[gg] = iter ;  
	}   
	//++++++++++++ END ITERATIONS person nn ++++++++++++++++++++++++    

	/////////////////////////////////////////////  
	// OUTPUT:  
	return Rcpp::List::create(  
			Rcpp::Named("theta") = theta ,  
			Rcpp::Named("setheta")=setheta ,  
			Rcpp::Named("Niter") = Niter  
		) ;  
}


//----------------------------------------------
     ////****** PRINT OUTPUT ****************************  
     // Rcpp::Rcout << "pij " << pij[0] << " " << pij[1] << " " <<   
     //    pij[2] << " " << pij[3] << " " << pij[4]        << std::endl ;  
     //Rcpp::Rcout << "ii  " << ii << " l1=" << l1 << " l2=" << l2 <<  
     //    "  lam1=" << lam1 << std::endl ;
	 //	Rcpp::Rcout << "nn ii  -- l2 " << l2 << std::endl ;  
//----------------------------------------------



///********************************************************************
///** mle_pcm_C
// [[Rcpp::export]]     
Rcpp::List mle_pcm_C( Rcpp::NumericMatrix dat, Rcpp::NumericMatrix dat_resp, Rcpp::NumericMatrix b, 
	Rcpp::NumericVector a, Rcpp::NumericVector maxK, Rcpp::NumericVector theta0, double conv, int maxiter )
{

	int N = dat.nrow() ;  
	int I = dat.ncol() ;  
	double eps=1e-10;
	double l1=0;  
	double l2=0;  
	double lam1=0;  
	double lam2=0;  
	double delta=0;  
	double sigdelta=0;  
	int K = b.ncol() + 1;  

	Rcpp::NumericVector theta(N) ;  
	Rcpp::NumericVector setheta(N) ;  
	Rcpp::NumericVector u(I) ;  
	Rcpp::NumericVector u_resp(I) ;  
	Rcpp::NumericVector pij0(K) ;  
	Rcpp::NumericVector pij(K);  
	Rcpp::NumericVector pij2(K);  
	Rcpp::NumericVector pij3(K) ;  
	Rcpp::NumericVector Niter(N) ;  
       
     double th1 = 0;  
     double t1=1;  
     double absdelta = 1000 ;  
     int iter = 0 ;  

     for (int nn=0; nn<N;nn++){  
         iter=0;  
         absdelta=100;  
         u = dat.row(nn) ;  
         u_resp = dat_resp.row(nn) ;  
       
         // iterations  
         th1 = theta0[nn] ;  
         iter=0;  
       
         //++++++++++++ BEGIN ITERATIONS person nn ++++++++++++++++++++++++    
         while( ( iter < maxiter ) & ( absdelta > conv)  ){  
             l1=0;  
             l2=0;          
             //************** begin ITEM ii **************  
             for (int ii=0; ii<I ; ii++){  
                 // calculate probabilities          
                 pij3[0] = 1;  
                 pij[0] = 0 ;  
                 t1=1;	                
                 for (int cc= 1 ; cc < maxK[ii] ; cc++){  
                     pij[cc] = pij[cc-1] + ( a[ii]* th1 - b(ii,cc-1) ) ;  
                     pij3[cc] = exp( pij[cc] ) ;  
                     t1 = pij3[cc] + t1 ;   
                             }  
                 for (int cc= 0 ; cc < maxK[ii] ; cc++){  
                     pij[cc] = pij3[cc] / t1 ;  
                             }  
                 // expected value  
                 lam1=0;  
                 lam2=0;  
                 for (int cc= 1;cc < maxK[ii];cc++){  
                     lam1 = lam1 + cc * pij[cc] ;  
                     lam2 = lam2 + cc*cc * pij[cc] ;  
                         }  
                 l1 = l1 + u_resp[ii] * a[ii] * ( u[ii] - lam1 );  
                 l2 = l2 - u_resp[ii] * a[ii]*a[ii] * ( lam2 - lam1 * lam1 ) ;    
                           
                         }  
         //************** end ITEM ii **************              
         // increment                  
         delta = - l1 / ( l2 + eps ) ;                                              
         sigdelta = 1 ;  // signum delta  
         if ( delta < 0 ){ sigdelta = -1 ; }  
         while( delta * sigdelta > 2){   delta = delta / 2 ; }   
         th1 = th1 + delta ;  
         absdelta = delta * sigdelta ;  
         iter ++ ; // increment iter                          
                   }  
                     
         theta[nn] = th1 ;            
         setheta[nn] = sqrt( - 1 / l2 ) ;  
         Niter[nn] = iter ;  
               }          
     //++++++++++++ END ITERATIONS person nn ++++++++++++++++++++++++    
                   
	/////////////////////////////////////////////  
	// OUTPUT:  
	return Rcpp::List::create(  
			Rcpp::Named("theta") = theta ,  
			Rcpp::Named("setheta")=setheta ,  
			Rcpp::Named("Niter") = Niter  
		) ;  
}



