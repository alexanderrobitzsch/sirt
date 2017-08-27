//// File Name: gooijercsntableaux.h
//// File Version: 1.08
//// File Last Change: 2014-03-12 18:03:32



// includes from the plugin
#include <RcppArmadillo.h>
#include <Rcpp.h>


using namespace Rcpp;




 
//**********************************************************************
// compute means and covariances of estimators obtained by Jackknife
Rcpp::NumericVector gta( Rcpp::NumericMatrix dat ){  
	
	int N=dat.nrow() ;
	int I=dat.ncol() ; // number of items
	int VV= I*(I-1)/ 2;
	int K=I+1 ;
	
	Rcpp::NumericVector score(N) ;
	Rcpp::NumericVector nk(I+1);
	Rcpp::NumericMatrix xik(I+1,I) ;
	Rcpp::NumericMatrix rijk(I+1,VV) ;
	Rcpp::NumericVector rijkmax(I+1) ;
	
	Rcpp::NumericVector stat(1) ;
	
	//****************
	// compute score distribution
	for( int nn=0;nn<N;nn++){
	for (int ii=0;ii<I;ii++){
		score[nn] += dat(nn,ii) ;
				}
		nk[ score[nn] ] ++ ;
			}
	nk = nk + 1E-20 ;
	
	//********
	// compute scorewise item means and counts
	
	for (int ii=0 ;ii<I;ii++){
	for (int nn=0 ;nn<N;nn++){
	    xik( score[nn] , ii ) += dat(nn,ii) ;	
		}  // end nn
	  for (int kk=0;kk<I+1;kk++){
	    xik(kk,ii ) = xik(kk,ii)/nk[kk] ;
		} // end kk
	}   // end ii
	
	//****************
	// compute conditional covariance
	
	int zz=0;
	// include something here!!!!!
	for (int ii =0;ii<I+1;ii++){
		rijkmax[ii] = rijkmax[ii] - 100 ;
				}
	
	for (int ii=0;ii<I;ii++){
	for (int jj=ii+1;jj<I;jj++){
	for (int nn=0;nn<N;nn++){
	    rijk( score[nn] , zz ) += dat(nn,ii)*dat(nn,jj) ;
				} // end nn
	for (int kk=0;kk<K;kk++){
	    rijk( kk , zz ) = rijk( kk , zz ) - nk[kk]*xik(kk,ii)*xik(kk,jj) ; 
	    rijk( kk , zz ) = rijk( kk , zz ) / nk[kk] ;
	    if ( rijk(kk,zz) > rijkmax[kk] ){
		    rijkmax[kk] = rijk(kk,zz) ;
				}
			}  // end kk
	    zz ++ ;
		} // end jj
	}  // end ii
	
	for (int kk=0;kk<K;kk++){
	   stat[0] += nk[kk] / N * rijkmax[kk] ;
			}
	
	return wrap( stat ) ;  
}

//***************************************************
// permutation of elements in dataset
Rcpp::NumericMatrix gooijer_permutation(Rcpp::NumericMatrix sampleM , int NS, int N ,
	Rcpp::NumericMatrix score_index ){
	Rcpp::NumericVector i1(2);
	int I = sampleM.ncol() ;
	int SS=score_index.nrow() ;
	int t1=0;
	int NS1=0;
	int n0=0;
	int n1=0;
//*** score_index table
//      ind          
//       1   1  1  1
// 64     2   5  4  3
// 1530   6   8  3  2
// 270    9  17  9  6
// 378   18  30 13  9
// 43    31  62 32 22
	
for ( int ss=0 ; ss < SS ; ss++){
	NS1 = score_index(ss,3 ) ;
	n0 = score_index(ss,0) - 1 ;
	n1 = score_index(ss,1) - 1 ;

// Rcpp::Rcout << "ss=" << ss  << std::endl ;	
	for (int ii =0 ; ii<I;ii++){ // variable ii
	for ( int zz=0;zz<NS1;zz++){
		i1 = Rcpp::floor( runif(2)*(n1-n0)+n0 ) ;
		t1 = sampleM( i1[0] ,ii ) ;
		sampleM( i1[0] ,ii)  = sampleM( i1[1] ,ii) ;
		sampleM( i1[1] , ii ) = t1 ;
			}  // end zz
		}  // end ii
	} // end ss
		
		
	return wrap( sampleM ) ;
		}
