//// File Name: sirt_rcpp_invariance_alignment.cpp
//// File Version: 2.575


// [[Rcpp::depends(RcppArmadillo)]]


// includes from the plugin
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[RcppNOinterfaces(r, cpp)]]

using namespace Rcpp;


///********************************************************************
///** sirt_rcpp_invariance_alignment_lambda_transformed
// [[Rcpp::export]]
Rcpp::NumericMatrix sirt_rcpp_invariance_alignment_lambda_transformed(
        Rcpp::NumericMatrix lambda, Rcpp::NumericVector psi0)
{
    int I = lambda.ncol();
    int G = lambda.nrow();
    Rcpp::NumericMatrix lambda1(G,I);
    for (int ii=0; ii<I; ii++){
        for (int gg=0; gg<G; gg++){
            lambda1(gg,ii) = lambda(gg,ii) / psi0[gg];
        }
    }
    //--- output
    return lambda1;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_invariance_alignment_nu_transformed
// [[Rcpp::export]]
Rcpp::NumericMatrix sirt_rcpp_invariance_alignment_nu_transformed(
        Rcpp::NumericMatrix nu, Rcpp::NumericMatrix lambda,
        Rcpp::NumericVector alpha0, Rcpp::NumericVector psi0,
        bool reparam)
{
    int I = nu.ncol();
    int G = nu.nrow();
    Rcpp::NumericMatrix nu1(G,I);
    for (int ii=0; ii<I; ii++){
        for (int gg=0; gg<G; gg++){
            if (reparam){
                nu1(gg,ii) = nu(gg,ii) - lambda(gg,ii) * alpha0[gg];
            } else {
                nu1(gg,ii) = nu(gg,ii) - lambda(gg,ii) / psi0[gg] * alpha0[gg];
            }
        }
    }
    //--- output
    return nu1;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_invariance_alignment_simplicity_function_value
// [[Rcpp::export]]
double sirt_rcpp_invariance_alignment_simplicity_function_value(
        Rcpp::CharacterVector type, double parm1, double parm2, double scale, double power,
        double eps)
{
    double diff_par = parm1 - parm2;
    double fval=0;
    double parm = std::pow(diff_par / scale, 2.0);
    if (power > 0){
        fval = std::pow( parm + eps, power );
    } else {
        double gamma0 = 50.0;
        fval = 2 / ( 1 + std::exp(-gamma0*std::sqrt( parm + eps)) ) - 1;
    }
    //--- output
    return fval;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_invariance_alignment_opt_fct
// [[Rcpp::export]]
Rcpp::List sirt_rcpp_invariance_alignment_opt_fct(
        Rcpp::NumericMatrix nu, Rcpp::NumericMatrix lambda,
        Rcpp::NumericVector alpha0, Rcpp::NumericVector psi0,
        Rcpp::IntegerMatrix group_combis, Rcpp::NumericMatrix wgt,
        Rcpp::NumericVector align_scale, Rcpp::NumericVector align_pow,
        double eps, Rcpp::NumericMatrix wgt_combi, Rcpp::CharacterVector type,
        bool reparam)
{
    int I = lambda.ncol();
    Rcpp::NumericMatrix lambda1 = sirt_rcpp_invariance_alignment_lambda_transformed(
                                        lambda, psi0);
    Rcpp::NumericMatrix nu1 = sirt_rcpp_invariance_alignment_nu_transformed(
                                    nu, lambda, alpha0, psi0, reparam);
    int GC = group_combis.nrow();
    double fopt=0;
    double parm1, parm2, parm3, parm4, fval;
    int gg1, gg2;
    for (int ii=0; ii<I; ii++){
        for (int cc=0; cc<GC; cc++){
            gg1 = group_combis(cc,0);
            gg2 = group_combis(cc,1);
            // simplicity function lambda
            parm1 = lambda1( gg1, ii);
            parm2 = lambda1( gg2, ii);
            fval = sirt_rcpp_invariance_alignment_simplicity_function_value(
                            type, parm1, parm2, align_scale[0], align_pow[0], eps);
            fopt += wgt_combi(cc,ii) * fval;
            // simplicity function nu
            parm3 = nu1( gg1, ii);
            parm4 = nu1( gg2, ii);
            fval = sirt_rcpp_invariance_alignment_simplicity_function_value(
                            type, parm3, parm4, align_scale[1], align_pow[1], eps);
            fopt += wgt_combi(cc,ii) * fval;
        }
    }
    //--- output
    return Rcpp::List::create(
            Rcpp::Named("fopt") = fopt,
            Rcpp::Named("lambda1") = lambda1,
            Rcpp::Named("nu1") = nu1
        );
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_invariance_alignment_simplicity_function_gradient
// [[Rcpp::export]]
double sirt_rcpp_invariance_alignment_simplicity_function_gradient(
        Rcpp::CharacterVector type, double parm1, double parm2, double scale, double power,
        double eps)
{
    double diff_par = parm1 - parm2;
    double fval=0;
    if (power > 0){
        fval = power * std::pow( std::pow(diff_par / scale, 2.0) + eps, power - 1.0 );
        fval = fval*2*diff_par/(scale*scale);
    } else {
        double gamma0 = 50.0;
        double g1 = ( 1 + std::exp(-gamma0*std::abs(diff_par/scale)) );
        double g2 = std::pow( std::pow(diff_par / scale, 2.0) + eps, -0.5 );
        g2 = g2*diff_par/(scale*scale);
        fval = 2*gamma0/g1 * ( 1 - 1 /g1 )* g2;
    }
    //--- output
    return fval;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_invariance_alignment_opt_grad
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_invariance_alignment_opt_grad(
        Rcpp::NumericMatrix nu, Rcpp::NumericMatrix lambda,
        Rcpp::NumericVector alpha0, Rcpp::NumericVector psi0,
        Rcpp::IntegerMatrix group_combis, Rcpp::NumericMatrix wgt,
        Rcpp::NumericVector align_scale, Rcpp::NumericVector align_pow,
        double eps, Rcpp::NumericMatrix wgt_combi, Rcpp::CharacterVector type,
        bool reparam)
{
    int I = lambda.ncol();
    int G = lambda.nrow();
    Rcpp::NumericMatrix lambda1 = sirt_rcpp_invariance_alignment_lambda_transformed(
                                        lambda, psi0);
    Rcpp::NumericMatrix nu1 = sirt_rcpp_invariance_alignment_nu_transformed(
                                    nu, lambda, alpha0, psi0, reparam);
    Rcpp::NumericVector grad(2*G);
    grad.fill(0);
    int GC = group_combis.nrow();
    double parm1, parm2, parm3, parm4, grad_val, temp1;
    int gg1, gg2, hh1, hh2;

    for (int ii=0; ii<I; ii++){
        for (int cc=0; cc<GC; cc++){
            gg1=group_combis(cc,0);
            gg2=group_combis(cc,1);
            hh1=gg1+G;
            hh2=gg2+G;
            //* simplicity function lambda
            parm1 = lambda1( gg1, ii);
            parm2 = lambda1( gg2, ii);
            grad_val = sirt_rcpp_invariance_alignment_simplicity_function_gradient(
                            type, parm1, parm2, align_scale[0], align_pow[0], eps);
            // parm=lambda1(gg,ii) = lambda(gg,ii) / psi0(gg,ii)
            // d parm / d psi = - lam / psi0^2 = - lam1 / psi0
            temp1 = wgt_combi(cc,ii) * grad_val;
            grad[hh1] -= temp1 * parm1 / psi0[gg1];
            grad[hh2] += temp1 * parm2 / psi0[gg2];

            //* simplicity function nu
            parm3 = nu1( gg1, ii);
            parm4 = nu1( gg2, ii);
            grad_val = sirt_rcpp_invariance_alignment_simplicity_function_gradient(
                            type, parm3, parm4, align_scale[1], align_pow[1], eps);
            // parm = nu1 = nu - lam / psi0 * alpha0
            temp1 = wgt_combi(cc,ii) * grad_val;
            grad[hh1] += temp1 * parm1 / psi0[gg1] * alpha0[gg1];
            grad[hh2] -= temp1 * parm2 / psi0[gg2] * alpha0[gg2];
            grad[gg1] -= temp1 * parm1;
            grad[gg2] += temp1 * parm2;
        }
    }
    //--- output
    return grad;
}
///********************************************************************


