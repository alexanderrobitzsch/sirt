//// File Name: sirt_rcpp_mgsem_functions.cpp
//// File Version: 0.365

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


///********************************************************************
///** computes the quadratic form t(y) %*% A %*% y
///** sirt_rcpp_mgsem_quadform
// [[Rcpp::export]]
double sirt_rcpp_mgsem_quadform(Rcpp::NumericMatrix y, Rcpp::NumericMatrix A)
{
    double val=0;
    int p=y.nrow();
    for (int ii=0; ii<p; ii++){
        val+=y[ii]*A(ii,ii)*y[ii];
    }
    for (int ii=0; ii<p-1; ii++){
        for (int jj=ii+1; jj<p; jj++){
            val+=2*y[ii]*A(ii,jj)*y[jj];
        }
    }
    //--- output
    return val;
}
///********************************************************************

///********************************************************************
///** computes the quadratic form t(y) %*% A %*% y with logical matrix
///** sirt_rcpp_mgsem_quadform_logical
// [[Rcpp::export]]
double sirt_rcpp_mgsem_quadform_logical(Rcpp::NumericMatrix y, Rcpp::NumericMatrix A,
            Rcpp::LogicalMatrix A_logical)
{
    double val=0;
    int p=y.nrow();
    for (int ii=0; ii<p; ii++){
        if (A_logical(ii,ii)){
            val+=y[ii]*A(ii,ii)*y[ii];
        }
    }
    for (int ii=0; ii<p-1; ii++){
        for (int jj=ii+1; jj<p; jj++){
            if (A_logical(ii,jj)){
                val+=2*y[ii]*A(ii,jj)*y[jj];
            }
        }
    }
    //--- output
    return val;
}
///********************************************************************

///********************************************************************
///** computes the trace of  A %*% B
///** this is sum(diag(A %*% B)) in pure R
///** sirt_rcpp_mgsem_trace_product
// [[Rcpp::export]]
double sirt_rcpp_mgsem_trace_product(Rcpp::NumericMatrix A, Rcpp::NumericMatrix B)
{
    double val=0;
    int p=A.nrow();
    for (int ii=0; ii<p; ii++){
        for (int jj=0; jj<p; jj++){
            val+=A(ii,jj)*B(jj,ii);
        }
    }
    //--- output
    return val;
}
///********************************************************************


///********************************************************************
///** computes the trace of  A %*% B using logical matrix for B entries
///** sirt_rcpp_mgsem_trace_product_logical
// [[Rcpp::export]]
double sirt_rcpp_mgsem_trace_product_logical(Rcpp::NumericMatrix A, Rcpp::NumericMatrix B,
            Rcpp::LogicalMatrix B_logical)
{
    double val=0;
    int p=A.nrow();
    for (int ii=0; ii<p; ii++){
        for (int jj=0; jj<p; jj++){
            if (B_logical(jj,ii)){
                val+=A(ii,jj)*B(jj,ii);
            }
        }
    }
    //--- output
    return val;
}
///********************************************************************

///********************************************************************
///** computes the trace of  A %*% B using logical matrix for B entries
///** sirt_rcpp_mgsem_loglike_derivative_parameters
// [[Rcpp::export]]
double sirt_rcpp_mgsem_loglike_derivative_parameters(
            Rcpp::NumericMatrix S1, Rcpp::NumericMatrix S3,
            Rcpp::NumericVector y, Rcpp::NumericMatrix Sigma_der,
            Rcpp::LogicalMatrix Sigma_der_logical)
{
    double val=0;
    int p=S1.nrow();
    for (int ii=0; ii<p; ii++){
        for (int jj=0; jj<p; jj++){
            // t2 <- -( t(y) %*% Sigma_der %*% y )[1,1]
            if (Sigma_der_logical(jj,ii)){
                val+=(S1(ii,jj)-S3(ii,jj))*Sigma_der(jj,ii);
                if (ii<jj){
                    val+=-2*Sigma_der(jj,ii)*y[ii]*y[jj];
                }
                if (ii==jj){
                    val+=-Sigma_der(jj,ii)*y[ii]*y[jj];
                }
            }
        }
    }
    //--- output
    return val;
}
///********************************************************************

///********************************************************************
///** computes model-implied covariance matrix
///** sirt_rcpp_mgsem_compute_cov
// [[Rcpp::export]]
Rcpp::NumericMatrix sirt_rcpp_mgsem_compute_cov(Rcpp::NumericMatrix LAM,
    Rcpp::NumericMatrix PHI, Rcpp::NumericMatrix PSI)
{
    int I=LAM.nrow();
    int D=LAM.ncol();
    Rcpp::NumericMatrix covmat(I,I);
    for (int ii=0; ii<I; ii++){
        for (int jj=ii; jj<I; jj++){
            covmat(ii,jj)=PSI(ii,jj);
            for (int dd=0; dd<D; dd++){
                for (int ee=0; ee<D; ee++){
                    covmat(ii,jj) += LAM(ii,dd)*PHI(dd,ee)*LAM(jj,ee);
                } //  end ee
            }    // end dd
            if (ii<jj){
                covmat(jj,ii)=covmat(ii,jj);
            }
        }    // end jj
    }    // end ii
    //--- output
    return covmat;
}
///********************************************************************


///********************************************************************
///** computes the quadratic form t(x) %*% y
///** sirt_rcpp_mgsem_sumproduct_logical
// [[Rcpp::export]]
double sirt_rcpp_mgsem_sumproduct_logical(Rcpp::NumericVector x, Rcpp::NumericVector y,
        Rcpp::LogicalVector y_logical)
{
    double val=0;
    int N=x.size();
    for (int ii=0; ii<N; ii++){
        if (y_logical[ii]){
            val += x[ii]*y[ii];
        }
    }
    //--- output
    return val;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_mgsem_vech_numeric
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_mgsem_vech_numeric(Rcpp::NumericMatrix A)
{
    int p=A.nrow();
    int n=p*(p+1)/2;
    Rcpp::NumericVector AV(n);
    int hh=0;
    for (int ii=0; ii<p; ii++){
        for (int jj=ii; jj<p;jj++){
            AV[hh]=A(ii,jj);
            hh++;
        }
    }
    //--- output
    return AV;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_mgsem_vech_logical
// [[Rcpp::export]]
Rcpp::LogicalVector sirt_rcpp_mgsem_vech_logical(Rcpp::LogicalMatrix A)
{
    int p=A.nrow();
    int n=p*(p+1)/2;
    Rcpp::LogicalVector AV(n);
    int hh=0;
    for (int ii=0; ii<p; ii++){
        for (int jj=ii; jj<p;jj++){
            AV[hh]=A(ii,jj);
            hh++;
        }
    }
    //--- output
    return AV;
}
///********************************************************************


///********************************************************************
///** sirt_rcpp_mgsem_eval_pen_lp_lasso
// [[Rcpp::export]]
double sirt_rcpp_mgsem_eval_pen_lp_lasso(double y, double eps_approx,
            double p, double regul_fac)
{
    double val=regul_fac*std::pow( y*y+eps_approx, p/2.0 );
    //--- output
    return val;
}
///********************************************************************


///********************************************************************
///** sirt_rcpp_mgsem_eval_pen_lp_scad
// [[Rcpp::export]]
double sirt_rcpp_mgsem_eval_pen_lp_scad(double y, double eps_approx,
            double p, double regul_fac, double a_scad)
{
    double x=std::pow( y*y+eps_approx, p/2.0 );
    double a=a_scad;
    double val=0;
    double lambda=regul_fac;

    // res <- ifelse( x < lambda, lambda * x, 0)
    if (x<lambda){
        val = lambda*x;
    }

    // res <- res + ifelse( ( x >=lambda ) & ( x < a*lambda),
    //                     - ( x^2 - 2*a*lambda*x+lambda^2) / ( 2*(a-1)),0 )
    if ( (x>=lambda) & (x<a*lambda) ){
        val = -( x*x - 2.0*a*lambda*x+lambda*lambda) / ( 2.0*(a-1));
    }

    // res <- res + ifelse( x>=a*lambda, (a+1)*lambda^2 / 2, 0 )
    if ( x>=a*lambda ){
        val = (a+1)*lambda*lambda / 2.0;
    }


    //--- output
    return val;
}
///********************************************************************


///********************************************************************
///** sirt_rcpp_mgsem_eval_pen_lp_lasso_deriv
// [[Rcpp::export]]
double sirt_rcpp_mgsem_eval_pen_lp_lasso_deriv(double y, double eps_approx,
            double p, double regul_fac)
{
    double val=regul_fac*std::pow( y*y+eps_approx, p/2.0-1.0 )*y;
    //--- output
    return val;
}
///********************************************************************


///********************************************************************
///** sirt_rcpp_mgsem_eval_lp_penalty
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_mgsem_eval_lp_penalty(Rcpp::NumericVector x,
        Rcpp::NumericVector fac, Rcpp::NumericVector n,
        double p, double eps_approx, bool deriv,
        Rcpp::CharacterVector pen_type, double a_scad, double h)
{

    int IP=x.size();
    Rcpp::NumericVector z(IP);
    double y1=0;
    double y2=0;

    for (int pp=0; pp<IP; pp++){
        if (fac[pp]>0){
            if (deriv){  // derivative
                if (pen_type[0]=="lasso"){
                    z[pp] = sirt_rcpp_mgsem_eval_pen_lp_lasso_deriv(x[pp], eps_approx,
                                p, fac[pp] );
                }
                if (pen_type[0]=="scad"){
                    y1 = sirt_rcpp_mgsem_eval_pen_lp_scad(x[pp]+h, eps_approx,
                                p, fac[pp], a_scad );
                    y2 = sirt_rcpp_mgsem_eval_pen_lp_scad(x[pp]-h, eps_approx,
                                p, fac[pp], a_scad );
                    z[pp] = (y1-y2)/(2.0*h);
                }

            } else {   // no derivative
                if (pen_type[0]=="lasso"){
                    z[pp] = sirt_rcpp_mgsem_eval_pen_lp_lasso(x[pp], eps_approx,
                                p, fac[pp] );
                }
                if (pen_type[0]=="scad"){
                    z[pp] = sirt_rcpp_mgsem_eval_pen_lp_scad(x[pp], eps_approx,
                                p, fac[pp], a_scad );
                }

            }  // end if deriv
            z[pp] = n[pp]*z[pp];
        }  // end fac
    }  // end pp

    //--- output
    return z;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_mgsem_eval_lpdiff_penalty
// [[Rcpp::export]]
double sirt_rcpp_mgsem_eval_lpdiff_penalty(Rcpp::NumericVector x,
            Rcpp::NumericMatrix fac, Rcpp::LogicalMatrix fac_logical,
            double p, double eps_approx, double a_scad,
            Rcpp::CharacterVector pen_type, Rcpp::NumericMatrix n)
{

    double z=0;
    double val=0;
    int IP=x.size();
    double y=0;

    for (int ii=0; ii<IP-1; ii++){
        for (int jj=ii+1; jj<IP; jj++){
            if (fac_logical(ii,jj)){
                y=std::abs(x[ii]-x[jj]);
                if (pen_type[0]=="lasso"){
                    // res <- (x^2+eps)^(p/2)
                    // val=fac(ii,jj)*std::pow( y*y+eps_approx, p/2.0 );
                    val = sirt_rcpp_mgsem_eval_pen_lp_lasso(y, eps_approx, p, fac(ii,jj) );
                }
                if (pen_type[0]=="scad"){
                    val = sirt_rcpp_mgsem_eval_pen_lp_scad(y, eps_approx, p, fac(ii,jj),
                                a_scad);
                }
                val=n(ii,jj)*val;
                z+=2*val;
            }  // end
        }  // end jj
    }  // end ii


    //--- output
    return z;
}
///********************************************************************


///********************************************************************
///** sirt_rcpp_mgsem_eval_lpdiff_penalty_deriv
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_mgsem_eval_lpdiff_penalty_deriv(Rcpp::NumericVector x,
        Rcpp::NumericMatrix fac, Rcpp::NumericMatrix n, Rcpp::LogicalMatrix fac_logical,
        double p, double eps_approx, double h, double a_scad,
        Rcpp::CharacterVector pen_type)
{
    int IP=x.size();
    Rcpp::NumericVector grad(IP);
    double ya=0;
    double yb=0;
    double y1=0;
    double y2=0;
    double n2=0;
    double h2=2.0*h;

    for (int pp=0; pp<IP; pp++){
        y1=0;
        y2=0;
        for (int jj=0; jj<IP; jj++){
            if (fac_logical(pp,jj)){
                ya=std::abs(x[pp]+h-x[jj]);
                yb=std::abs(x[pp]-h-x[jj]);
                n2=2.0*n(pp,jj);
                if (pen_type[0]=="lasso"){
                    y1+=n2*sirt_rcpp_mgsem_eval_pen_lp_lasso(ya, eps_approx, p,
                                            fac(pp,jj) );
                    y2+=n2*sirt_rcpp_mgsem_eval_pen_lp_lasso(yb, eps_approx, p,
                                            fac(pp,jj) );
                }
                if (pen_type[0]=="scad"){
                    y1+=n2*sirt_rcpp_mgsem_eval_pen_lp_scad(ya, eps_approx, p,
                                        fac(pp,jj),    a_scad);
                    y2+=n2*sirt_rcpp_mgsem_eval_pen_lp_scad(yb, eps_approx, p,
                                        fac(pp,jj), a_scad );
                }
            }  // end
        }  // end jj
        grad[pp] = (y1-y2)/h2;
    }  // end pp

    //--- output
    return grad;
}
///********************************************************************
