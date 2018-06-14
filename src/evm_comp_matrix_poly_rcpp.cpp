//// File Name: evm_comp_matrix_poly_rcpp.cpp
//// File Version: 3.617

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;




///////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List firsteigenvalsirt2(arma::mat X, int maxit, double conv, double K)
{
    double lambda_temp ;
    double lambda =0 ;
    double lambdadiff = 1000;
    double lambda_old ;
    Rcpp::NumericVector lambda_est(2);
    int iter = 0 ;

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



//**********************************************************************
// compute means and covariances of estimators obtained by Jackknife
// [[Rcpp::export]]
Rcpp::List parameters_jackknife( Rcpp::NumericMatrix PARS )
{
    int VV=PARS.nrow() ;
    int JJ=PARS.ncol() ;
    Rcpp::NumericVector PARS_means(VV);
    Rcpp::NumericMatrix PARS_vcov(VV,VV);
    Rcpp::NumericMatrix PARS_(VV,JJ);
    for (int jj=0; jj<JJ; jj++){
        PARS_(_,jj) = PARS(_,jj);
    }

    double tmp3=0;
    // compute row means
    for (int vv=0;vv<VV;vv++){
        tmp3=0;
        for (int jj=0;jj<JJ;jj++){
            tmp3+= PARS_(vv,jj) ;
        }
        PARS_means[vv] = tmp3 / JJ ;
        PARS_(vv,_) = PARS_(vv,_) - PARS_means[vv] ;
    }
    // compute covariance
    for (int vv1=0;vv1<VV;vv1++){
        for (int vv2=vv1;vv2<VV;vv2++){
            for (int jj=0;jj<JJ;jj++){
                PARS_vcov(vv1,vv2) += PARS_(vv1,jj)*PARS_(vv2,jj) ;
            } // end jj
            PARS_vcov(vv1,vv2) = PARS_vcov(vv1,vv2) * (JJ-1) / JJ ;
            if (vv1!=vv2){
                PARS_vcov(vv2,vv1) = PARS_vcov(vv1,vv2);
            }
        }
    }
    //------- output
    return Rcpp::List::create(
                Rcpp::Named("PARS_means")= PARS_means,
                Rcpp::Named("PARS_vcov")= PARS_vcov
            );
}



//************************************************************
// Eigenvector method
// [[Rcpp::export]]
Rcpp::List evm_aux( arma::mat B , int I , int powD ,
    int maxit, double conv, double K )
{
    int n1 = B.n_rows;
    int n2 = B.n_cols;
    arma::mat TMP = B( arma::span(0, n1-1), arma::span(0, n2-1) );

    arma::mat D = arma::zeros(I,I);
    for ( int hh=0;hh<(powD-1) ; hh++){
        TMP = TMP * B ;
    }
    //  calculate ratios in D
    for (int ii=0;ii<I;ii++){
        D(ii,ii) = 1 ;
        for (int jj=0;jj<I;jj++){
            if ( ii!= jj){
                D(ii,jj) = TMP(jj,ii) / TMP(ii,jj);
            }
        }
    }
    // compute first eigenvalue
    Rcpp::List res11 = firsteigenvalsirt2( D, maxit, conv, K) ;
    Rcpp::NumericVector u = res11["u"] ;
    double tmp1= 0 ;
    // center difficulty vector
    for ( int ii = 0 ; ii<I ; ii++){
        u[ii] = log( u[ii] ) ;
        tmp1 += u[ii] ;
    }
    Rcpp::NumericVector b = u - tmp1 / I ;

    // compute eigenvalue and consistency index
    Rcpp::NumericVector tmp2 = res11["lambda1"] ;
    double lambda = tmp2[0] ;
    double cons_index = ( lambda - I )   / ( I -1 ) ;

    return Rcpp::List::create(
                Rcpp::Named("lambda")= lambda ,
                Rcpp::Named("D") = D ,
                Rcpp::Named("cons_index") = cons_index ,
                Rcpp::Named("b") = b
            ) ;
}


// Rcpp::Rcout << "a100  " <<   std::endl ;

//***************************************
// Choppin's row averaging approach
// [[Rcpp::export]]
Rcpp::NumericVector choppin_rowaveraging( arma::mat B , int I , double priorweight)
{
    Rcpp::NumericVector b_ra(I) ;
    arma::mat TMP2= arma::zeros(I,I);

    int n1 = B.n_rows;
    int n2 = B.n_cols;
    arma::mat B_ = B( arma::span(0, n1-1), arma::span(0, n2-1) );
    B_ = B_ + priorweight ;

    //  calculate ratios in D
    for (int ii=0;ii<I;ii++){
        for (int jj=0;jj<I;jj++){
            if (ii!=jj){
                TMP2(ii,jj) = log(  B_(jj,ii)  /  B_(ii,jj)  )  ;
            }
        }
    }
    for (int ii=0;ii<I;ii++){
        for (int jj=0;jj<I;jj++){
            b_ra[ii] += TMP2(ii,jj ) ;
        }
        b_ra[ii] = b_ra[ii] / I ;
    }
    //--- output
    return  b_ra;
}
//**********************************************



///********************************************************************
///** evm_comp_matrix_poly
// [[Rcpp::export]]
Rcpp::List evm_comp_matrix_poly( Rcpp::NumericMatrix dat ,
    Rcpp::NumericMatrix dat_resp , Rcpp::NumericVector weights ,
    int JJ , Rcpp::NumericVector jackunits, int powD, int progress_,
    Rcpp::NumericMatrix row_index, Rcpp::NumericMatrix col_index )
{
    int N = dat.nrow() ;
    int I = row_index.nrow() ;
    double K = I ;
    double conv= 0.0001;
    int maxit=100;

    // matrix for dichotomous responses
    arma::mat B = arma::zeros( I , I  );
    arma::mat Bjack = arma::zeros( I , I*JJ ) ;
    arma::mat B2 = arma::zeros( I , I ) ;

    Rcpp::NumericMatrix b_ra_jack(I , JJ ) ;
    Rcpp::NumericMatrix b_evm_jack(I , JJ ) ;
    Rcpp::NumericVector lambda_jack(JJ) ;
    Rcpp::NumericVector cons_index_jack(JJ) ;
    Rcpp::NumericVector lambda_jack_na(JJ) ;

    //*****************************************
    // start counting pairwise comparisons

    if ( progress_==1){
        Rcpp::Rcout << "*** Create pairwise comparison matrix " << std::flush <<  std::endl ;
    }

    int ii=0;
    int jj=0;
    for (int rr=0;rr<I;rr++){
        for (int cc=0;cc<I;cc++){
            ii = row_index(rr,0) ;
            jj = col_index(cc,0) ;
            if ( ii != jj ){
                for (int nn=0;nn<N;nn++){
                    if ( dat_resp(nn,ii) * dat_resp(nn,jj) == 1 ){
                        if ( ( dat(nn,ii)==row_index(rr,1) ) & ( ( dat(nn,jj)==col_index(cc,1))) ){
                            B(rr,cc) += weights[nn];
                            Bjack(rr,cc + I*jackunits[nn] ) += weights[nn] ;
                        }
                    }  // end if dat_resp(nn,ii)==1 & dat_resp(nn,jj)==1
                }    // end nn
            } // end if ii != jj
        } // end cc
    } // end rr

    ///***************************************
    // eigenvector method

    if ( progress_==1){
        Rcpp::Rcout << "*** Perform eigenvector method "  << std::flush <<  std::endl ;
    }

    //***  original data
    Rcpp::List res0 = evm_aux(B , I , powD , maxit, conv, K);
    Rcpp::NumericVector lambda = res0["lambda"];
    Rcpp::NumericMatrix D = res0["D"];
    Rcpp::NumericVector cons_index = res0["cons_index"] ;
    Rcpp::NumericVector b = res0["b"];
    Rcpp::NumericVector tmp1, tmp2;

    //*** jackknife
    int JJadj = JJ ;
    for (int jj=0;jj<JJ;jj++){
        B2=0*B2 ;
        B2 = B - Bjack.submat( arma::span(0,I-1), arma::span(jj*I, I-1+jj*I) );
        if (JJ==1){
            B2 = B;
        }
        Rcpp::List res2 = evm_aux(B2 , I , powD , maxit, conv, K) ;
        Rcpp::NumericVector v1 = res2["b"] ;
        b_evm_jack(_,jj) = v1 ;
        tmp1 = res2["lambda"];
        lambda_jack[jj] = tmp1[0] ;
        if ( R_IsNaN( tmp1[0] ) ){
            lambda_jack_na[jj] = 1 ;
            JJadj = JJadj - 1 ;
        }
        tmp2 = res2["cons_index"] ;
        cons_index_jack[jj] = tmp2[0] ;
    }

    //************************************************
    // statistical inference
    // create matrix with all parameters
    int VV = B.n_rows + 2 ;
    Rcpp::NumericMatrix PARS( VV , JJadj ) ;
    int jj2=0;
    for (int jj=0;jj<JJ;jj++){
        if ( lambda_jack_na[jj] == 0 ){
            for (int pp=0;pp<VV-2;pp++){
                PARS(pp,jj2) = b_evm_jack(pp,jj) ;
                PARS(VV-2,jj2 ) = lambda_jack[jj] ;
                PARS(VV-1,jj2 ) = cons_index_jack[jj] ;
            }
            jj2 ++;
        }
    }
    // inference based on Jackknife
    Rcpp::List res1 =  parameters_jackknife(PARS);
    Rcpp::NumericMatrix PARS_vcov = res1["PARS_vcov"];
    Rcpp::NumericVector PARS_means = res1["PARS_means"];

    //*************************************************
    // OUTPUT
    return Rcpp::List::create(
                Rcpp::Named("D") = D ,
                Rcpp::Named("B") = B ,
                Rcpp::Named("b_evm") = b ,
                Rcpp::Named("b_evm_jack") = b_evm_jack ,
                Rcpp::Named("lambda") = lambda ,
                Rcpp::Named("lambda_jack") = lambda_jack ,
                Rcpp::Named("lambda_jack_na") = lambda_jack_na ,
                Rcpp::Named("cons_index") = cons_index,
                Rcpp::Named("cons_index_jack") = cons_index_jack,
                Rcpp::Named("JJ") = JJ,
                Rcpp::Named("JJadj") = JJadj ,
                Rcpp::Named("PARS_means") = PARS_means ,
                Rcpp::Named("PARS_vcov") = PARS_vcov
            ) ;
}

