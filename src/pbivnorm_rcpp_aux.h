//// File Name: pbivnorm_rcpp_aux.h
//// File Version: 0.52


// includes from the plugin
#include <RcppArmadillo.h>
#include <Rcpp.h>


using namespace Rcpp;



 
//**********************************************************************
// bivariate normal distribution


Rcpp::NumericVector pbivnorm2_rcpp( Rcpp::NumericVector x_ , Rcpp::NumericVector y_  , 
	Rcpp::NumericVector rho_  ) 
{
	
	Rcpp::NumericVector x(x_) ;
	Rcpp::NumericVector y(y_) ;
	Rcpp::NumericVector rho1(rho_) ;	
	
	double a=x[0] ;
	double b=y[0] ;
	double rho=rho1[0] ;	
	//    a <- x
	//    b <- y
	//    a11 <- a1 <- - a
	double a11 = -a ;
	double a1 = -a ;
	//    b1 <- - b
	double b1=-b;
	
	//	# APPROXIMATION for negative correlations
	//	ind.neg <- which( rho < 0 )
	//	if ( length(ind.neg) > 0 ){
	//		rho[ind.neg] <- - rho[ind.neg]
	//		b1[ind.neg] <- -b1[ind.neg]
	//					}
	int ind_neg = 0;
	if ( rho < 0 ){
		ind_neg = 1 ;
		rho = -rho ;
		b1 = -b1 ;
			}	
			
	//	# APPROX. (ii)
	//	ind2 <- which( a1 < 0 & b1 < 0 )
	//	if ( length(ind2) > 0 ){
	//		a1[ ind2 ] <- - a1[ind2] 
	//		b1[ ind2 ] <- - b1[ind2] 
	//				}
	int ind2 = 0 ;
	if ( (a1 < 0) & ( b1 < 0 ) ){
		ind2 = 1;
		a1 = -a1 ;
		b1 = -b1 ;
			}
	//	# APPROX. (i):
	//	# a > 0 and b > 0 => a > b
	//	ind <- which( a1 < b1 )	
	//	t1 <- b1
	//	if ( length(ind) > 0 ){ 
	//		b1[ ind ] <- a1[ind]
	//		a1[ind] <- t1[ind]
	//			}
	// int ind=0;
	double t1 = b1 ;
	if ( a1 < b1 ){
	//	ind = 1 ;
		b1 = a1 ;
		a1 = t1 ;
			}
	//    t1 <- pnorm( - a1 )
	t1 = Rf_pnorm5( - a1 , 0.0, 1.0, 1, 0); 
	//    mu <- dnorm( a1 ) / pnorm( - a1 )
	double mu = Rf_dnorm4(a1,0.0,1.0,FALSE) / t1 ;
	//    xi <-  ( rho * mu - b1 ) / sqrt( 1 - rho^2 ) 
	double rho2 = rho*rho ;
	double xi = ( rho * mu - b1 ) / sqrt( 1 - rho2) ;
	//    sig2 <- 1 + a1*mu - mu^2
	double sig2 = 1+a1*mu - mu*mu ;
	// prob1 <- t1 * ( pnorm(  xi ) - 1/2 * rho^2 / ( 1 - rho^2 ) * xi * dnorm( xi ) * sig2  )
	double prob1 = t1 * (  Rf_pnorm5(  xi , 0.0, 1.0, 1, 0) - 0.5*rho2 / ( 1 - rho2 ) *
		xi * Rf_dnorm4(xi,0.0,1.0,FALSE) * sig2 ) ;	
	//	if ( length(ind2) > 0 ){
	//		prob1[ind2] <- 1 - pnorm( -a1[ind2] ) - pnorm( -b1[ind2] ) + prob1[ind2] 
	//			}
	if ( ind2 == 1 ){
	    prob1 = 1 - Rf_pnorm5(  -a1 , 0.0, 1.0, 1, 0) - Rf_pnorm5(  -b1 , 0.0, 1.0, 1, 0) +
			prob1 ;
			}
	//	if ( length(ind.neg) > 0 ){
	//		prob1[ind.neg] <- pnorm( - a11[ind.neg] ) - prob1[ ind.neg]
	//					}	
	if ( ind_neg == 1){
	   prob1 = Rf_pnorm5(  -a11 , 0.0, 1.0, 1, 0) - prob1 ;
			}
	Rcpp::NumericVector prob2(1);
	prob2[0] = prob1 ;
	return wrap(prob2) ;

}


//**********************************************************************
// bivariate normal density


Rcpp::NumericVector dmvnorm_2dim_rcpp( Rcpp::NumericVector x_, 
   	   Rcpp::NumericVector y_, Rcpp::NumericVector rho_){
        
     //  x , y , rho   
     Rcpp::NumericVector x1_(x_) ;  
     Rcpp::NumericVector y1_(y_) ;     
     Rcpp::NumericVector rho1_(rho_) ;          
     double x = x1_[0] ;
     double y = y1_[0] ;
     double rho = rho1_[0] ;         
     double tmp1 = 0 ;  
     double pi1 =  3.14159265358979323846 ;    
     tmp1 = x*x - 2*rho*x*y+y*y ;  
     tmp1 = tmp1 / ( 2*(1-rho*rho) ) ;  
     tmp1 = exp( - tmp1 ) ;   
     tmp1 = tmp1 / ( 2*pi1*sqrt(1-rho*rho) ) ;          
     Rcpp::NumericVector d1(1);   
     d1[0] = tmp1 ;         
     return wrap(d1) ;         
}


//**********************************************************************
// calculation estimating equation polychoric correlation


Rcpp::NumericVector polychoric2_estequation( Rcpp::NumericMatrix frtab ,
	int maxK , Rcpp::NumericVector rho , Rcpp::NumericVector thresh1n ,
	Rcpp::NumericVector thresh2n , int maxK1 , int maxK2 
  	   ) 

{
// INPUT:
// frtab , thresh1n , thresh2n , rho , maxK1 , maxK2             
   Rcpp::NumericMatrix pi_ij(maxK+2,maxK+2 );
   Rcpp::NumericMatrix li_ij(maxK+2,maxK+2 );
   Rcpp::NumericMatrix phip_ij(maxK+2,maxK+2 );
   Rcpp::NumericMatrix phid_ij(maxK+2,maxK+2 );
   Rcpp::NumericVector tmp1ii(1) ;
   Rcpp::NumericVector tmp1jj(1) ;
   Rcpp::NumericVector tmp2(1);
   double eps2 = 1e-15 ;
	// compute distribution and density
	for ( int ii=0 ; ii < maxK1+1 ; ii++){
	for ( int jj=0 ; jj < maxK2+1 ; jj++){	
		tmp1ii[0] = thresh1n[ii] ;
		tmp1jj[0] = thresh2n[jj] ;	
		phip_ij(ii,jj) = pbivnorm2_rcpp(  tmp1ii , tmp1jj , rho )[0] ; 
		phid_ij(ii,jj) = dmvnorm_2dim_rcpp( tmp1ii , tmp1jj , rho )[0] ;
					}
				}	

	// compute pi_ij and estimating equation
	Rcpp::NumericVector llest(1) ;
	for (int ii=1 ; ii <maxK1+1;ii++){
	for (int jj=1 ; jj < maxK2+1 ; jj++){
	   pi_ij(ii,jj) = phip_ij(ii,jj) - phip_ij(ii-1,jj) - phip_ij(ii,jj-1) + phip_ij(ii-1,jj-1);
	   li_ij(ii,jj) = phid_ij(ii,jj) - phid_ij(ii-1,jj) - phid_ij(ii,jj-1) + phid_ij(ii-1,jj-1);
	   li_ij(ii,jj) =  frtab(ii-1,jj-1) * li_ij(ii,jj) / ( pi_ij(ii,jj) + eps2 ) ;	  
	   llest[0] += li_ij(ii,jj) ;	   
					}
				}
     return wrap(llest) ;
     
     	   }

     	   

//**********************************************************************
// calculation polychoric correlation for one item pair



Rcpp::List polychoric2_itempair( Rcpp::NumericVector v1 , Rcpp::NumericVector v2 ,
	int maxK , int maxiter
  	   ) 

{

maxK = maxK + 1 ;
int CC = v1.size() ;
double eps=1e-5;
double numdiffparm =1e-6 ;
double conv=1e-10;

//////////*********** FREQUENCIES ************************//////

// calculate frequency table
int Ntotal = 0 ; // Ntotal 
Rcpp::NumericMatrix frtab(maxK,maxK);
for (int cc = 0 ; cc < CC ; cc++){
	if ( ( ! R_IsNA( v1[cc] ) ) & ( ! R_IsNA( v2[cc] ) ) ){
		frtab( v1[cc] , v2[cc] ) ++ ;
		Ntotal ++ ;		
		}
	}
// set values to .5 if they are zero
for (int cc1=0;cc1 < maxK;cc1++){
for (int cc2=0;cc2 < maxK;cc2++){
if (frtab(cc1,cc2) < 1){
	frtab(cc1,cc2) = .5 ;
			}
}
}
	
// compute marginals from frequency table ;
Rcpp::NumericVector thresh1(maxK+2) ;
Rcpp::NumericVector thresh1n(maxK+2) ;
Rcpp::NumericVector thresh2(maxK+2) ;
Rcpp::NumericVector thresh2n(maxK+2) ;
int maxK1=0;
int maxK2=0;
double tmp1 = 0 ;


thresh1n[0] = eps ;
thresh2n[0] = eps ;
thresh1n[maxK+1] = 1-eps ;
thresh2n[maxK+1] = 1-eps ;



for ( int zz = 0 ; zz < maxK ; zz++){
   tmp1=0;
	for (int cc=0;cc<maxK;cc++){
		tmp1 += frtab(zz,cc) ;
			}
	if (tmp1 > 0 ){ maxK1 = zz ; }
	thresh1[zz+1]=thresh1[zz]+tmp1 ;
	thresh1n[zz+1] = thresh1[zz+1] / ( Ntotal + eps ) ;
	if ( thresh1n[zz+1] > 1 - 1E-10 ){ thresh1n[zz+1] = 1 - eps ; }
	if ( thresh1n[zz+1] == 0 ){ thresh1n[zz+1] = eps ; }	
}

for ( int zz = 0 ; zz < maxK ; zz++){
   tmp1=0;
	for (int cc=0;cc<maxK;cc++){
		tmp1 += frtab(cc,zz) ;
			}
	if (tmp1 > 0 ){ maxK2 = zz ; }
	thresh2[zz+1]=thresh2[zz]+tmp1 ;
	thresh2n[zz+1] = thresh2[zz+1] / ( Ntotal+eps ) ;
	if ( thresh2n[zz+1] > 1 - 1E-10 ){ thresh2n[zz+1] = 1 - eps ; }
	if ( thresh2n[zz+1] == 0 ){ thresh2n[zz+1] = eps ; }	
}	

thresh1n = Rcpp::qnorm( thresh1n , 0.0 , 1.0 ) ;
thresh2n = Rcpp::qnorm( thresh2n , 0.0 , 1.0 ) ;

maxK1 = maxK1+1 ;
maxK2 = maxK2+1 ;

//////////*********** ALGORITHM POLYCHORIC ************************//////

Rcpp::NumericVector rho(1) ;  // initial correlation = 0 ;
Rcpp::NumericVector rho1(1) ;

double ll0 , ll1 ;
double incr , der ;
int iter =0 ;
double aincr=1;

while ( ( iter < maxiter ) & ( aincr > conv) & ( Ntotal > 0 )  ){

	if (rho[0] > 1 - numdiffparm ){
		rho[0] = rho[0] - numdiffparm*1.5 ;	
				}		
	rho1[0] = rho[0] + numdiffparm ;	
	ll0 = polychoric2_estequation( frtab ,maxK,rho,thresh1n,thresh2n,maxK1,maxK2)[0] ; 
	ll1 = polychoric2_estequation( frtab ,maxK,rho1,thresh1n,thresh2n,maxK1,maxK2)[0] ;
	der = (ll1-ll0)/numdiffparm ;
	incr = ll0 / der ;
//	Rcout << "within iter " << iter << " incr " << incr << std::endl ;	
	rho[0] = rho[0] - incr ;
	if (rho[0] > 1- eps ){ 
		rho[0] = .90 + .01 * iter ;		
			}
	if (rho[0] < -1+ eps ){ 
		rho[0] = -.90 - .01 * iter ;		
			}			
	if (rho[0] < - 1 + eps){ rho[0] = -1 + eps ; }
	if (rho[0] > 1 ){ 
		rho[0] = 1 - eps ;		
			}	
	aincr = incr ;
	if (incr < 0 ){ aincr = - incr ; }
 // Rcout << "iter " << iter << "rho " << rho0[0]  << std::endl;	
	iter ++ ;
		}

Rcpp::NumericVector maxK11(1);		
Rcpp::NumericVector maxK21(1);

maxK11[0] = maxK1-1 ;
maxK21[0] = maxK2-1 ;

Rcpp::NumericVector Ntot1(1);
Ntot1[0] = Ntotal ;

//*************************************************    
// OUTPUT            
            
return Rcpp::List::create(  
    Rcpp::_["thresh1"] = thresh1n , 
    Rcpp::_["thresh2"] = thresh2n ,    
    Rcpp::_["rho"] = rho  , 
    Rcpp::_["maxK1"] = maxK11 , 
    Rcpp::_["maxK2"] = maxK21 ,
    Rcpp::_["Ntotal"]=Ntot1 
    ) ;  
	
}
