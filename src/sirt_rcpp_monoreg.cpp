//// File Name: sirt_rcpp_monoreg.cpp
//// File Version: 1.05



// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


// user includes


///********************************************************************
///** sirt_rcpp_monoreg_rowwise
// [[Rcpp::export]]
Rcpp::NumericMatrix sirt_rcpp_monoreg_rowwise( Rcpp::NumericMatrix YM,
    Rcpp::NumericMatrix WM )
{
    //********************************************
    // Adapted code from fdrtool::monoreg function
    // contained in the fdrtool package
    // original author: Korbinian Strimmer
    //********************************************

    // define row and column numbers
    int NR=YM.nrow();
    int NC=YM.ncol();
    int nn=NC;

    // create output ghat matrix
    Rcpp::NumericMatrix ghat (NR,NC);
    Rcpp::NumericMatrix k (NR,NC);
    Rcpp::NumericMatrix gew (NR,NC);
    double neu;
    int zz, c, j;

    for (zz=0; zz<NR; zz++){ // begin for zz
        c=0;
        j=0;
        nn = NC;
        k(zz,c) = 0;
        gew(zz,c) = WM(zz,0);
        ghat(zz,c) = YM(zz,0);

        for (j=1; j<nn; j++){   // begin for j ...
            c = c+1;
            k(zz,c) = j;
            gew(zz,c) = WM(zz,j);
            ghat(zz,c) = YM(zz,j);
            /* c is at least 1 as nn is > 1 */
            while (ghat(zz,c-1) >= ghat(zz,c)){
                neu = gew(zz,c)+gew(zz,c-1);
                ghat(zz,c-1) = ghat(zz,c-1)+(gew(zz,c)/neu)*(ghat(zz,c)-ghat(zz,c-1) );
                gew(zz,c-1) = neu;
                c = c-1;
                if (c==0) break;
            }
        }     // end for j ...
        //##########
        while (nn >= 1){
            for (j=k(zz,c); j<nn; j++){
                ghat(zz,j) = ghat(zz,c);
            }
            nn = k(zz,c);
            c = c-1;
        }
    } // end for zz

    ///---- OUTPUT
    return ghat;
}
///********************************************************************
