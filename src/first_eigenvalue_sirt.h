//// File Name: first_eigenvalue_sirt.h
//// File Version: 2.22


// includes from the plugin
#include <RcppArmadillo.h>
#include <Rcpp.h>


using namespace Rcpp;


// declarations
// extern "C" {
// Rcpp::List firsteigenvalsirt(arma::mat X, int maxit, double conv, double K) ;
// }


Rcpp::List firsteigenvalsirt(arma::mat X, int maxit, double conv, double K)
{
	double lambda_temp;
	double lambda=0;
	double lambdadiff=1000;
	double lambda_old;
	Rcpp::NumericVector lambda_est(2);
	int iter = 0;

	//**********
	// set matrices
	arma::mat Xz ;
	arma::colvec z(K) ;
	double temp1 = 1 / sqrt( K ) ;
	for (int ii=0;ii<K;ii++){
		z[ii] = temp1 ;
	}

	///////////////////////////////
	/// algorithm

	while ( ( iter < maxit ) & ( lambdadiff > conv) ){
		lambda_old = lambda ;
		Xz = arma::mat( X * z ) ;
		lambda_temp = 0 ;
		for (int ii=0;ii<K;ii++){
			lambda_temp += Xz[ii]*Xz[ii] ;
		}
		lambda = sqrt( lambda_temp ) ;
		lambdadiff = lambda - lambda_old ;
		if ( lambdadiff < 0 ){ lambdadiff = - lambdadiff ; }
		z = Xz / lambda ;
		iter ++ ;
	}

	lambda_est[0] = lambda ;

	////////////////////////////////////
	// OUTPUT:
	return Rcpp::List::create(
				Rcpp::Named("u")=z, 
				Rcpp::Named("lambda1")=lambda_est 
			) ;
}


    // Rcpp::Rcout << "Iteration =" << iter << " lambda = " << lambda << std::endl ; 
    // Rcpp::Rcout << " difference = " << lambdadiff << std::endl ; 



