//// File Name: polychoric2_tetrachoric2_rcpp.cpp
//// File Version: 3.328


// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


//**********************************************************************
// bivariate normal distribution
///** pbivnorm2_rcpp
// [[Rcpp::export]]
Rcpp::NumericVector pbivnorm2_rcpp( Rcpp::NumericVector x_, Rcpp::NumericVector y_,
    Rcpp::NumericVector rho_  )
{
    Rcpp::NumericVector x(x_);
    Rcpp::NumericVector y(y_);
    Rcpp::NumericVector rho1(rho_);

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
    if ( (a1 < 0) & ( b1 < 0 ) ){
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
    double t1 = b1;
    if ( a1 < b1 ){
        b1 = a1;
        a1 = t1;
    }
    //    t1 <- pnorm( - a1 )
    t1 = ::Rf_pnorm5( -a1, 0.0, 1.0, 1, 0);
    //    mu <- dnorm( a1 ) / pnorm( - a1 )
    double mu = ::Rf_dnorm4(a1, 0.0, 1.0, FALSE) / t1;
    //    xi <-  ( rho * mu - b1 ) / sqrt( 1 - rho^2 )
    double rho2 = rho*rho;
    double xi = ( rho * mu - b1 ) / std::sqrt( 1 - rho2);
    //    sig2 <- 1 + a1*mu - mu^2
    double sig2 = 1 + a1*mu - mu*mu;
    // prob1 <- t1 * ( pnorm(  xi ) - 1/2 * rho^2 / ( 1 - rho^2 ) * xi * dnorm( xi ) * sig2  )
    double prob1 = t1 * ( ::Rf_pnorm5(  xi, 0.0, 1.0, 1, 0) - 0.5*rho2 / ( 1 - rho2 ) *
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


const double pi1 = 3.1415926535897;

//**********************************************************************
// bivariate normal density
///** dmvnorm_2dim_rcpp
// [[Rcpp::export]]
Rcpp::NumericVector dmvnorm_2dim_rcpp( Rcpp::NumericVector x_,
        Rcpp::NumericVector y_, Rcpp::NumericVector rho_)
{
    //  x, y, rho
    Rcpp::NumericVector x1_(x_);
    Rcpp::NumericVector y1_(y_);
    Rcpp::NumericVector rho1_(rho_);
    double x = x1_[0];
    double y = y1_[0];
    double rho = rho1_[0];
    double tmp1 = 0;
    tmp1 = x*x - 2*rho*x*y + y*y;
    tmp1 = tmp1 / ( 2*(1-rho*rho) );
    tmp1 = std::exp( - tmp1 );
    tmp1 = tmp1 / ( 2*pi1*std::sqrt(1-rho*rho) );
    Rcpp::NumericVector d1(1);
    d1[0] = tmp1;
    return d1;
}
//**********************************************************************

//**********************************************************************
// calculation estimating equation polychoric correlation
///** polychoric2_estequation
// [[Rcpp::export]]
Rcpp::NumericVector polychoric2_estequation( Rcpp::NumericMatrix frtab,
    int maxK, Rcpp::NumericVector rho, Rcpp::NumericVector thresh1n,
    Rcpp::NumericVector thresh2n, int maxK1, int maxK2 )
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
            phip_ij(ii,jj) = pbivnorm2_rcpp( tmp1ii, tmp1jj, rho )[0];
            phid_ij(ii,jj) = dmvnorm_2dim_rcpp( tmp1ii, tmp1jj, rho )[0];
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
///** polychoric2_itempair
// [[Rcpp::export]]
Rcpp::List polychoric2_itempair( Rcpp::NumericVector v1, Rcpp::NumericVector v2,
    int maxK_, int maxiter)
{
    int maxK = maxK_ + 1;
    int CC = v1.size();
    double eps=1e-5;
    double h =1e-6;
    double conv=1e-10;

    //////////*********** FREQUENCIES ************************//////
    // calculate frequency table
    int Ntotal = 0; // Ntotal
    Rcpp::NumericMatrix frtab(maxK,maxK);
    for (int cc = 0; cc < CC; cc++){
        if ( ( ! R_IsNA( v1[cc] ) ) & ( ! R_IsNA( v2[cc] ) ) ){
            frtab( v1[cc], v2[cc] ) ++;
            Ntotal ++;
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
    Rcpp::NumericVector rho(1);  // initial correlation = 0;
    Rcpp::NumericVector rho1(1);

    double ll0, ll1;
    double incr, der;
    int iter =0;
    double aincr=1;
    while ( ( iter < maxiter ) & ( aincr > conv) & ( Ntotal > 0 )  ){
        if (rho[0] > 1 - h ){
            rho[0] = rho[0] - h*1.5;
        }
        rho1[0] = rho[0] + h;
        ll0 = polychoric2_estequation( frtab,maxK,rho,thresh1n,thresh2n,maxK1,maxK2)[0];
        ll1 = polychoric2_estequation( frtab,maxK,rho1,thresh1n,thresh2n,maxK1,maxK2)[0];
        der = (ll1-ll0)/h;
        incr = ll0 / der;
        rho[0] = rho[0] - incr;
        if (rho[0] > 1- eps ){
            rho[0] = .90 + .01 * iter;
        }
        if (rho[0] < -1+ eps ){
            rho[0] = -.90 - .01 * iter;
        }
        if (rho[0] < -1 + eps){
            rho[0] = -1 + eps;
        }
        if (rho[0] > 1 ){
            rho[0] = 1 - eps;
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
    //*************************************************
    // OUTPUT
    return Rcpp::List::create(
            Rcpp::Named("thresh1") = thresh1n,
            Rcpp::Named("thresh2") = thresh2n,
            Rcpp::Named("rho") = rho,
            Rcpp::Named("maxK1") = maxK11,
            Rcpp::Named("maxK2") = maxK21,
            Rcpp::Named("Ntotal") = Ntot1
        );
}


// Rcout << "iter " << iter << "rho " << rho0[0]  << std::endl;

///********************************************************************
///** tetrachoric2_rcpp_aux
// [[Rcpp::export]]
Rcpp::NumericVector tetrachoric2_rcpp_aux( Rcpp::NumericMatrix dfr,
    double h, int maxiter )
{
    int ZZ = dfr.nrow();
    Rcpp::NumericVector x(1);
    Rcpp::NumericVector y(1);
    Rcpp::NumericVector rho(1);
    Rcpp::NumericVector p11(1);
    Rcpp::NumericVector rhores(ZZ);
    Rcpp::NumericVector L0(1), L0h(1);
    double L1;
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
        while ( ( iter < maxiter ) & ( aincr > maxincr ) ){
            L0 = pbivnorm2_rcpp( x, y, rho );
            L0h = pbivnorm2_rcpp( x, y, rho+h );
            L1 = ( L0h[0] - L0[0] ) / h;
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


    // Rcout << "prob1 " << prob1 << std::endl;


///********************************************************************
///** polychoric2_aux_rcpp
// [[Rcpp::export]]
Rcpp::List polychoric2_aux_rcpp( Rcpp::NumericMatrix dat,
    int maxK, int maxiter )
{
    int I = dat.ncol();
    int maxK3 = maxK + 3;
    Rcpp::NumericVector v1, tmp1, tmp2;
    Rcpp::NumericVector v2;
    Rcpp::List res;
    Rcpp::NumericMatrix rho(I,I);
    Rcpp::NumericMatrix Nobs(I,I);
    Rcpp::NumericMatrix thresh(I,maxK3);
    Rcpp::NumericVector maxcat(I);
    Rcpp::NumericVector tmp3;
    Rcpp::NumericVector Ntot_used(I);

    for (int ii=0; ii < I-1;ii++){
        for (int jj=ii+1;jj<I;jj++){
            v1 = dat(_,ii);
            v2 = dat(_,jj);
            res=polychoric2_itempair( v1, v2, maxK, maxiter );
            tmp1 = res["rho"];
            rho(ii,jj) = tmp1[0];
            rho(jj,ii) = rho(ii,jj);
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

    //*************************************************
    // OUTPUT
    return Rcpp::List::create(
            Rcpp::Named("rho") = rho,
            Rcpp::Named("thresh")=thresh,
            Rcpp::Named("maxcat") = maxcat,
            Rcpp::Named("Nobs")=Nobs,
            Rcpp::Named("Ntot_used") = Ntot_used
        );
}

