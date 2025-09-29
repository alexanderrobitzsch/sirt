//// File Name: sirt_rcpp_polychoric2.cpp
//// File Version: 3.454


// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[RcppNOinterfaces(r, cpp)]]

using namespace Rcpp;

// [include_header_file]
#include "sirt_rcpp_pbvnorm.h"

#include <pbv.h>
// # include "c:/Users/sunpn563/Dropbox/Eigene_Projekte/R-Routinen/IRT-Functions/pbv_package/0.3/pbv_work/inst/include/pbv.h"

const double pi1 = 3.1415926535897;

//**********************************************************************
// bivariate normal distribution
///** sirt_rcpp_pbivnorm2
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_pbivnorm2( Rcpp::NumericVector x, Rcpp::NumericVector y,
    Rcpp::NumericVector rho1  )
{
    double a=x[0];
    double b=y[0];
    double rho=rho1[0];
    //    a <- x
    //    b <- y
    //    a11 <- a1 <- - a
    double a11 = -a;
    double a1 = -a;
    //    b1 <- - b
    double b1=-b;

    //    # APPROXIMATION for negative correlations
    //    ind.neg <- which( rho < 0 )
    //    if ( length(ind.neg) > 0 ){
    //        rho[ind.neg] <- - rho[ind.neg]
    //        b1[ind.neg] <- -b1[ind.neg]
    //    }
    int ind_neg = 0;
    if ( rho < 0 ){
        ind_neg = 1;
        rho = -rho;
        b1 = -b1;
    }
    //    # APPROX. (ii)
    //    ind2 <- which( a1 < 0 & b1 < 0 )
    //    if ( length(ind2) > 0 ){
    //        a1[ ind2 ] <- - a1[ind2]
    //        b1[ ind2 ] <- - b1[ind2]
    //                }
    int ind2 = 0;
    if ( (a1 < 0) && ( b1 < 0 ) ){
        ind2 = 1;
        a1 = -a1;
        b1 = -b1;
    }
    //    # APPROX. (i):
    //    # a > 0 and b > 0 => a > b
    //    ind <- which( a1 < b1 )
    //    t1 <- b1
    //    if ( length(ind) > 0 ){
    //        b1[ ind ] <- a1[ind]
    //        a1[ind] <- t1[ind]
    //            }
    // int ind=0;
    double eps1 = 1e-3;
    double t1 = b1;
    if ( a1 < b1 ){
        b1 = a1;
        a1 = t1;
    }
    //    t1 <- pnorm( - a1 )
    t1 = ::Rf_pnorm5( -a1, 0.0, 1.0, 1, 0);
    if (t1 < eps1){ t1 = eps1; }
    //    mu <- dnorm( a1 ) / pnorm( - a1 )
    double mu = ::Rf_dnorm4(a1, 0.0, 1.0, FALSE) / t1;
    //    xi <-  ( rho * mu - b1 ) / sqrt( 1 - rho^2 )
    double rho2 = rho*rho;
    double rho21 = 1 + 1e-5 - rho2;
    double xi = ( rho * mu - b1 ) / std::sqrt(rho21);
    //    sig2 <- 1 + a1*mu - mu^2
    double sig2 = 1 + a1*mu - mu*mu;
    // prob1 <- t1 * ( pnorm(  xi ) - 1/2 * rho^2 / ( 1 - rho^2 ) * xi * dnorm( xi ) * sig2  )
    double prob1 = t1 * ( ::Rf_pnorm5(  xi, 0.0, 1.0, 1, 0) - 0.5*rho2 / rho21  *
                        xi * ::Rf_dnorm4(xi, 0.0, 1.0, FALSE) * sig2 );
    //    if ( length(ind2) > 0 ){
    //        prob1[ind2] <- 1 - pnorm( -a1[ind2] ) - pnorm( -b1[ind2] ) + prob1[ind2]
    //            }
    if ( ind2 == 1 ){
        prob1 = 1 - ::Rf_pnorm5( -a1, 0.0, 1.0, 1, 0) - ::Rf_pnorm5( -b1, 0.0, 1.0, 1, 0) + prob1;
    }
    //    if ( length(ind.neg) > 0 ){
    //        prob1[ind.neg] <- pnorm( - a11[ind.neg] ) - prob1[ ind.neg]
    //    }
    if ( ind_neg == 1){
        prob1 = ::Rf_pnorm5( -a11, 0.0, 1.0, 1, 0) - prob1;
    }
    Rcpp::NumericVector prob2(1);
    prob2[0] = prob1;
    return prob2;
}
//**********************************************************************


//**********************************************************************
// bivariate normal density
///** sirt_rcpp_dmvnorm_2dim
// [[Rcpp::export]]
double sirt_rcpp_dmvnorm_2dim( double x, double y, double rho)
{
    double tmp1 = 0;
    double eps = 1e-10;
    tmp1 = x*x - 2*rho*x*y + y*y;
    tmp1 = tmp1 / ( 2*(1+eps-rho*rho) );
    tmp1 = std::exp( - tmp1 );
    tmp1 = tmp1 / ( 2*pi1*std::sqrt(1+eps-rho*rho) );
    return tmp1;
}
//**********************************************************************


//**********************************************************************
///** sirt_rcpp_polychoric2_pbivnorm
// [[Rcpp::export]]
double sirt_rcpp_polychoric2_pbivnorm( Rcpp::NumericVector x,
        Rcpp::NumericVector y, Rcpp::NumericVector rho, int use_pbv)
{
    double val=0;
    if (use_pbv==2){
        val = pbv::pbv_rcpp_pbvnorm0( x[0], y[0], rho[0]);
    }
    if (use_pbv==1){
        val = sirt_rcpp_pbvnorm0( x[0], y[0], rho[0]);
    }
    if (use_pbv==0){
        val = sirt_rcpp_pbivnorm2( x, y, rho )[0];
    }
    return val;
}
//**********************************************************************



//**********************************************************************
// calculation estimating equation polychoric correlation
///** sirt_rcpp_polychoric2_estimating_equation
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_polychoric2_estimating_equation(
    Rcpp::NumericMatrix frtab, int maxK, Rcpp::NumericVector rho,
    Rcpp::NumericVector thresh1n, Rcpp::NumericVector thresh2n, int maxK1,
    int maxK2, int use_pbv )
{
    // INPUT:
    // frtab, thresh1n, thresh2n, rho, maxK1, maxK2
    Rcpp::NumericMatrix pi_ij(maxK+2,maxK+2);
    Rcpp::NumericMatrix li_ij(maxK+2,maxK+2);
    Rcpp::NumericMatrix phip_ij(maxK+2,maxK+2);
    Rcpp::NumericMatrix phid_ij(maxK+2,maxK+2);
    Rcpp::NumericVector tmp1ii(1);
    Rcpp::NumericVector tmp1jj(1);
    Rcpp::NumericVector tmp2(1);
    double eps2 = 1e-15;
    // compute distribution and density
    for ( int ii=0; ii < maxK1+1; ii++){
        for ( int jj=0; jj < maxK2+1; jj++){
            tmp1ii[0] = thresh1n[ii];
            tmp1jj[0] = thresh2n[jj];
            phip_ij(ii,jj) = sirt_rcpp_polychoric2_pbivnorm( tmp1ii, tmp1jj, rho, use_pbv );
            phid_ij(ii,jj) = sirt_rcpp_dmvnorm_2dim( tmp1ii[0], tmp1jj[0], rho[0] );
        }
    }
    // compute pi_ij and estimating equation
    Rcpp::NumericVector llest(1);
    for (int ii=1; ii < maxK1+1;ii++){
        for (int jj=1; jj < maxK2+1; jj++){
            pi_ij(ii,jj) = phip_ij(ii,jj) - phip_ij(ii-1,jj) - phip_ij(ii,jj-1) + phip_ij(ii-1,jj-1);
            li_ij(ii,jj) = phid_ij(ii,jj) - phid_ij(ii-1,jj) - phid_ij(ii,jj-1) + phid_ij(ii-1,jj-1);
            li_ij(ii,jj) = frtab(ii-1,jj-1) * li_ij(ii,jj) / ( pi_ij(ii,jj) + eps2 );
            llest[0] += li_ij(ii,jj);
        }
    }
    return llest;
}
//**********************************************************************


//**********************************************************************
// calculation polychoric correlation for one item pair
///** sirt_rcpp_polychoric2_est_itempair
// [[Rcpp::export]]
Rcpp::List sirt_rcpp_polychoric2_est_itempair( Rcpp::NumericVector v1,
    Rcpp::NumericVector v2, int maxK_, int maxiter, int use_pbv,
    double conv, double rho_init, Rcpp::NumericVector weights )
{
    int maxK = maxK_ + 1;
    int CC = v1.size();
    double eps=1e-5;
    double h=1e-5;

    //////////*********** FREQUENCIES ************************//////
    // calculate frequency table
    int Ntotal = 0; // Ntotal
    Rcpp::NumericMatrix frtab(maxK,maxK);
    for (int cc = 0; cc < CC; cc++){
        if ( ( ! R_IsNA( v1[cc] ) ) && ( ! R_IsNA( v2[cc] ) ) ){
            frtab( v1[cc], v2[cc] ) += weights[cc];
            Ntotal += weights[cc];
        }
    }
    // set values to .5 if they are zero
    for (int cc1=0;cc1 < maxK;cc1++){
        for (int cc2=0;cc2 < maxK;cc2++){
            if (frtab(cc1,cc2) < 1){
                frtab(cc1,cc2) = .5;
            }
        }
    }

    // compute marginals from frequency table;
    Rcpp::NumericVector thresh1(maxK+2);
    Rcpp::NumericVector thresh1n(maxK+2);
    Rcpp::NumericVector thresh2(maxK+2);
    Rcpp::NumericVector thresh2n(maxK+2);
    int maxK1=0;
    int maxK2=0;
    double tmp1 = 0;
    thresh1n[0] = eps;
    thresh2n[0] = eps;
    thresh1n[maxK+1] = 1-eps;
    thresh2n[maxK+1] = 1-eps;

    for ( int zz = 0; zz < maxK; zz++){
        tmp1=0;
        for (int cc=0;cc<maxK;cc++){
            tmp1 += frtab(zz,cc);
        }
        if (tmp1 > 0 ){ maxK1 = zz; }
        thresh1[zz+1]=thresh1[zz]+tmp1;
        thresh1n[zz+1] = thresh1[zz+1] / ( Ntotal + eps );
        if ( thresh1n[zz+1] > 1 - conv ){ thresh1n[zz+1] = 1 - eps; }
        if ( thresh1n[zz+1] == 0 ){ thresh1n[zz+1] = eps; }
    }

    for ( int zz = 0; zz < maxK; zz++){
        tmp1=0;
        for (int cc=0;cc<maxK;cc++){
            tmp1 += frtab(cc,zz);
        }
        if (tmp1 > 0 ){ maxK2 = zz; }
        thresh2[zz+1] = thresh2[zz]+tmp1;
        thresh2n[zz+1] = thresh2[zz+1] / ( Ntotal+eps );
        if ( thresh2n[zz+1] > 1 - conv ){ thresh2n[zz+1] = 1 - eps; }
        if ( thresh2n[zz+1] == 0 ){ thresh2n[zz+1] = eps; }
    }

    thresh1n = Rcpp::qnorm( thresh1n, 0.0, 1.0 );
    thresh2n = Rcpp::qnorm( thresh2n, 0.0, 1.0 );
    maxK1 = maxK1+1;
    maxK2 = maxK2+1;

    //////////*********** ALGORITHM POLYCHORIC ************************//////
    Rcpp::NumericVector rho(1);
    Rcpp::NumericVector rho1(1);
    Rcpp::NumericVector rho2(1);
    rho1[0] = rho_init;
    rho[0] = rho_init;

    double ll0=0;
    double ll1=0;
    double incr=0;
    double der=0;
    int iter =0;
    double aincr=1;
    while ( ( iter < maxiter ) && ( aincr > conv) && ( Ntotal > 0 )  ){
        if (rho[0] > 1 - h ){
            rho[0] = rho[0] - h*1.5;
        }
        rho1[0] = rho[0] + h;
        rho2(0) = rho(0) - h;
        ll0 = sirt_rcpp_polychoric2_estimating_equation( frtab, maxK, rho2, thresh1n,
                            thresh2n, maxK1, maxK2, use_pbv)[0];
        ll1 = sirt_rcpp_polychoric2_estimating_equation( frtab, maxK, rho1, thresh1n,
                            thresh2n, maxK1, maxK2, use_pbv)[0];
        der = (ll1-ll0)/(2*h);
        incr = ll0 / der;
        rho[0] = rho[0] - incr;
        if (rho[0] > 1-eps ){
            rho[0] = .90 + .01*iter;
        }
        if (rho[0] < -1+eps ){
            rho[0] = -.90 - .01*iter;
        }
        if (rho[0] < -1+eps){
            rho[0] = -1+eps;
        }
        if (rho[0] > 1 ){
            rho[0] = 1-eps;
        }
        aincr = std::abs(incr);
        iter ++;
    }

    Rcpp::NumericVector maxK11(1);
    Rcpp::NumericVector maxK21(1);
    maxK11[0] = maxK1-1;
    maxK21[0] = maxK2-1;
    Rcpp::NumericVector Ntot1(1);
    Ntot1[0] = Ntotal;

    //--- OUTPUT
    return Rcpp::List::create(
            Rcpp::Named("thresh1") = thresh1n,
            Rcpp::Named("thresh2") = thresh2n,
            Rcpp::Named("rho") = rho,
            Rcpp::Named("maxK1") = maxK11,
            Rcpp::Named("maxK2") = maxK21,
            Rcpp::Named("Ntotal") = Ntot1,
            Rcpp::Named("iter") = iter,
            Rcpp::Named("der") = der
        );
}
///********************************************************************


// Rcout << "iter " << iter << "rho " << rho0[0]  << std::endl;
// Rcout << "iter " << iter << " incr=" << incr  << std::endl;

///********************************************************************
///** sirt_rcpp_polychoric2
// [[Rcpp::export]]
Rcpp::List sirt_rcpp_polychoric2( Rcpp::NumericMatrix dat, int maxK,
    int maxiter, int use_pbv, double conv, Rcpp::IntegerMatrix rho_init,
    Rcpp::NumericVector weights)
{
    int I = dat.ncol();
    int maxK3 = maxK + 3;
    Rcpp::NumericVector v1, v2, tmp1, tmp2, tmp3;
    Rcpp::List res;
    Rcpp::NumericMatrix rho(I,I);
    Rcpp::IntegerMatrix iter(I,I);
    Rcpp::NumericMatrix Nobs(I,I);
    Rcpp::NumericMatrix thresh(I,maxK3);
    Rcpp::NumericVector maxcat(I);
    Rcpp::NumericVector Ntot_used(I);

    for (int ii=0; ii < I-1;ii++){
        for (int jj=ii+1;jj<I;jj++){
            v1 = dat(_,ii);
            v2 = dat(_,jj);
            res = sirt_rcpp_polychoric2_est_itempair( v1, v2, maxK, maxiter,
                            use_pbv, conv, rho_init(jj,ii), weights );
            tmp1 = res["rho"];
            rho(ii,jj) = tmp1[0];
            rho(jj,ii) = rho(ii,jj);
            iter(jj,ii) = res["iter"];
            tmp1 = res["Ntotal"];
            Nobs(ii,jj) = tmp1[0];
            Nobs(jj,ii) = Nobs(ii,jj);
            Ntot_used(ii) += tmp1[0];
            Ntot_used(jj) += tmp1[0];
            tmp2 = res["thresh1"];
            for (int uu=0;uu<maxK3;uu++){
                thresh(ii,uu) += tmp2[uu]*Nobs(ii,jj);
            }
            tmp2 = res["thresh2"];
            for (int uu=0;uu<maxK3;uu++){
                thresh(jj,uu) += tmp2[uu]*Nobs(ii,jj);
            }
            tmp3 = res["maxK1"];
            if (tmp3[0] > maxcat[ii] ){ maxcat[ii] = tmp3[0]; }
            tmp3 = res["maxK2"];
            if (tmp3[0] > maxcat[jj] ){ maxcat[jj] = tmp3[0]; }
        }
    }

    for (int ii=0; ii < I;ii++){
        for (int uu=0;uu<maxK3; uu++){
            thresh(ii,uu) = thresh(ii,uu) / Ntot_used(ii);
            if ( uu > maxcat[ii]  ){
                thresh(ii,uu) = 99;
            }
        }
        rho(ii,ii) = 1;
    }

    //--- OUTPUT
    return Rcpp::List::create(
            Rcpp::Named("rho") = rho,
            Rcpp::Named("thresh")=thresh,
            Rcpp::Named("maxcat") = maxcat,
            Rcpp::Named("Nobs")=Nobs,
            Rcpp::Named("Ntot_used") = Ntot_used,
            Rcpp::Named("iter") = iter
        );
}
///********************************************************************



///********************************************************************
///** sirt_rcpp_tetrachoric2_rcpp
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_tetrachoric2_rcpp( Rcpp::NumericMatrix dfr,
    double h, int maxiter )
{
    int ZZ = dfr.nrow();
    Rcpp::NumericVector x(1);
    Rcpp::NumericVector y(1);
    Rcpp::NumericVector rho(1);
    Rcpp::NumericVector p11(1);
    Rcpp::NumericVector rhores(ZZ);
    Rcpp::NumericVector L0(1), L0h(1);
    double L1=0;
    double aincr=1;
    double incr=0;
    double maxincr = 1e-9;
    double eps=1e-10;
    int iter=0;

    for (int zz = 0; zz<ZZ;zz++){
        iter=0;
        aincr=1;
        x[0] = dfr(zz,10);
        y[0] = dfr(zz,11);
        rho[0] = dfr(zz,15);
        p11[0] = dfr(zz,7);
        while ( ( iter < maxiter ) && ( aincr > maxincr ) ){
            L0 = sirt_rcpp_pbivnorm2( x, y, rho-h );
            L0h = sirt_rcpp_pbivnorm2( x, y, rho+h );
            L1 = ( L0h[0] - L0[0] ) / (2*h);
            // Newton-Raphson step according Divgi (1979)
            incr = ( L0[0] - p11[0] ) / ( L1 + eps );
            rho[0] = rho[0] - incr;
            iter ++;
            aincr = std::abs(incr);
        }
        rhores[zz] = rho[0];
    }
    return rhores;
}
///********************************************************************



// Rcout << "prob1 " << prob1 << std::endl;
