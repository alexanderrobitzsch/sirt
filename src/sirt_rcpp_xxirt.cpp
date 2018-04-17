//// File Name: sirt_rcpp_xxirt.cpp
//// File Version: 0.22



// [[Rcpp::depends(RcppArmadillo)]]

// #include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


///********************************************************************
///** sirt_rcpp_xxirt_compute_posterior_expected_counts
// [[Rcpp::export]]           
Rcpp::NumericMatrix sirt_rcpp_xxirt_compute_posterior_expected_counts( 
		Rcpp::LogicalMatrix dat1_resp_gg, Rcpp::NumericMatrix p_aj_xi_gg )
{
	int N = dat1_resp_gg.nrow();
	int I = dat1_resp_gg.ncol();
	int TP = p_aj_xi_gg.ncol(); 	

	Rcpp::NumericMatrix nij(I, TP);        
	double val=0;

	//*** loop over items and categories     
	for (int ii=0;ii<I;ii++){
		for (int tt=0;tt<TP;tt++){
			val=0;
			for (int nn=0;nn<N;nn++){
				if ( dat1_resp_gg(nn,ii) ){
					val +=  p_aj_xi_gg(nn,tt) ;
				}
			}  // end nn
			nij(ii,tt) = val;
		}   // end tt
	}   // end ii  

	//--- OUTPUT              
	return nij ;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_xxirt_compute_likelihood
// [[Rcpp::export]]           
Rcpp::NumericMatrix sirt_rcpp_xxirt_compute_likelihood( 
		Rcpp::IntegerMatrix dat, Rcpp::LogicalMatrix dat_resp_bool,
		Rcpp::NumericMatrix probs, int TP, int maxK )
{
	int N = dat.nrow();
	int I = dat.ncol();
	Rcpp::NumericMatrix p_xi_aj(N, TP);

	for (int nn=0;nn<N;nn++){
		for (int tt=0;tt<TP;tt++){
			p_xi_aj(nn,tt) = 1 ;
		}
		for (int ii=0;ii<I;ii++){
			if ( dat_resp_bool(nn,ii) ){
				for (int tt=0;tt<TP;tt++){
					p_xi_aj(nn,tt) = p_xi_aj(nn,tt) * probs(ii , dat(nn,ii) + tt*maxK );
				}
			}
		} 
	}

	//---- OUTPUT
	return p_xi_aj ;    
}
///********************************************************************

// if ( ! R_IsNA( resp(nn,ii) ) ){
