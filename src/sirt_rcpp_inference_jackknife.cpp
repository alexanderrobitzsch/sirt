//// File Name: sirt_rcpp_inference_jackknife.cpp
//// File Version: 0.03



// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;




//**********************************************************************
// compute means and covariances of estimators obtained by Jackknife
//* sirt_rcpp_inference_jackknife
// [[Rcpp::export]]
Rcpp::List sirt_rcpp_inference_jackknife( Rcpp::NumericMatrix PARS )
{
    int VV=PARS.nrow();
    int JJ=PARS.ncol();
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
            tmp3+= PARS_(vv,jj);
        }
        PARS_means[vv] = tmp3 / JJ;
        PARS_(vv,_) = PARS_(vv,_) - PARS_means[vv];
    }
    // compute covariance
    for (int vv1=0;vv1<VV;vv1++){
        for (int vv2=vv1;vv2<VV;vv2++){
            for (int jj=0;jj<JJ;jj++){
                PARS_vcov(vv1,vv2) += PARS_(vv1,jj)*PARS_(vv2,jj);
            } // end jj
            PARS_vcov(vv1,vv2) = PARS_vcov(vv1,vv2) * (JJ-1) / JJ;
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
//**********************************************************************
