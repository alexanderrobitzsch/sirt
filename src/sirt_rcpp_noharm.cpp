//// File Name: sirt_rcpp_noharm.cpp
//// File Version: 3.454

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


///********************************************************************
///** sirt_rcpp_noharm_sirt_optim_fn_rcpp
// [[Rcpp::export]]
double sirt_rcpp_noharm_sirt_optim_fn_rcpp(Rcpp::NumericMatrix gamma_val,
    Rcpp::NumericVector delta, int I, Rcpp::NumericMatrix wgtm, Rcpp::NumericMatrix pm,
    Rcpp::NumericMatrix b0_jk, Rcpp::NumericMatrix b1_jk, Rcpp::NumericMatrix b2_jk,
    Rcpp::NumericMatrix b3_jk)
{
    double val=0;
    double x_ij, g2, g3, pm_exp;
    for (int ii=0; ii<I-1; ii++){
        for (int jj=ii+1; jj<I; jj++){
            if (wgtm(ii,jj) >0 ){
                x_ij = gamma_val(ii,jj) / std::sqrt( delta[ii] * delta[jj] );
                g2 = x_ij*x_ij;
                g3 = g2*x_ij;
                pm_exp = b0_jk(ii,jj) + b1_jk(ii,jj)*x_ij + b2_jk(ii,jj)*g2 + b3_jk(ii,jj)*g3;
                val = val + wgtm(ii,jj)*std::pow( pm(ii,jj) - pm_exp, 2);
            }
        }
    }
    //--- output
    return val;
}
///********************************************************************


///********************************************************************
///** sirt_rcpp_noharm_sirt_optim_gr_rcpp
// [[Rcpp::export]]
Rcpp::List sirt_rcpp_noharm_sirt_optim_gr_rcpp(Rcpp::NumericMatrix gamma_val,
    int NH, int I, Rcpp::NumericMatrix wgtm, Rcpp::NumericMatrix pm,
    Rcpp::NumericMatrix b0_jk, Rcpp::NumericMatrix b1_jk, Rcpp::NumericMatrix b2_jk,
    Rcpp::NumericMatrix b3_jk, int npar, Rcpp::IntegerVector pt_matid,
    Rcpp::IntegerVector pt_index, Rcpp::IntegerVector pt_row,
    Rcpp::IntegerVector pt_col, Rcpp::NumericMatrix FP, Rcpp::NumericMatrix Fmat,
    Rcpp::NumericMatrix Pmat, Rcpp::NumericMatrix Psimat)
{
    //---- derivative with respect to gamma diagonale
    Rcpp::NumericMatrix grad_gamma_diag(I, npar);
    Rcpp::LogicalMatrix grad_gamma_diag_bool(I, npar);

    int mat_hh, par_index_hh, row, col;
    double der;

    for (int ii=0; ii<I; ii++){
        for (int hh=0; hh<NH; hh++){
            //     mat_hh <- parm_table_free[hh,"mat"]
            //    par_index_hh <- parm_table_free[hh,"index"]
            //    row <- parm_table_free[hh,"row"]
            //    col <- parm_table_free[hh,"col"]
            mat_hh = pt_matid[hh];
            par_index_hh = pt_index[hh];
            row = pt_row[hh];
            col = pt_col[hh];
            // # F
            if (mat_hh==1){
                if (row==ii){
                    der = 2*FP(ii,col);
                    grad_gamma_diag(ii,par_index_hh) += der;
                    grad_gamma_diag_bool(ii,par_index_hh) = TRUE;
                }
            }
            // # P
            if (mat_hh==2){
                der = 0;
                if (row==col){
                    der = Fmat(ii,col)*Fmat(ii,col);
                    grad_gamma_diag(ii,par_index_hh) += der;
                } else {
                    der = Fmat(ii,row)*Fmat(ii,col);
                    grad_gamma_diag(ii,par_index_hh) += der;
                }
                if ( std::abs(der) > 0 ){
                    grad_gamma_diag_bool(ii,par_index_hh) = TRUE;
                }
            }
        } // end hh
    }  // end ii

    //----- transformations
    Rcpp::NumericVector delta(I);
    Rcpp::NumericMatrix grad_gamma_diag1(I,npar);
    for (int ii=0; ii<I; ii++){
        delta[ii] = 1 + gamma_val(ii,ii);
        for (int pp=0; pp<npar; pp++){
            grad_gamma_diag1(ii,pp) = grad_gamma_diag(ii,pp) / delta[ii];
        }
    }

    //--------------------------------
    // derivative with respect to item pairs
    double gii, gjj, gij, za, val, x_ij, b0, b1, b2, b3, g1, g2, g3;
    double pm_exp, pm_exp_der, temp1;

    Rcpp::NumericVector grad0(npar);
    Rcpp::NumericVector grad(npar);

    for (int ii=0; ii<I-1; ii++){
        for (int jj=ii+1; jj<I; jj++){
            if (wgtm(ii,jj) > 0 ){
                grad.fill(0);
                for (int hh=0; hh<NH; hh++){
                    mat_hh = pt_matid[hh];
                    par_index_hh = pt_index[hh];
                    row = pt_row[hh];
                    col = pt_col[hh];
                    // # F
                    if (mat_hh == 1){
                        der = 0;
                        if (row==ii){
                            der = FP(jj,col);
                            grad[par_index_hh] += der;
                        }
                        if (row==jj){
                            der = FP(ii,col);
                            grad[par_index_hh] += der;
                        }
                    }
                    // # P
                    if (mat_hh == 2){
                        if (row==col){
                            der = Fmat(ii,col)*Fmat(jj,col);
                            grad[par_index_hh] += der;
                        } else {
                            der = Fmat(ii,row)*Fmat(jj,col);
                            grad[par_index_hh] += der;
                        }
                    }
                    // # Psi
                    if (mat_hh == 3){
                        if (row==ii){
                            if (col==jj){
                                der = 1;
                                grad[par_index_hh] += der;
                            }
                        }
                    }
                } // end hh

                gii = (1+gamma_val(ii,ii));
                gjj = (1+gamma_val(jj,jj));
                gij = gamma_val(ii,jj);
                za = std::sqrt(gii*gjj);
                val = gij / za;
                // grad <- grad / za;
                // grad <- grad - 0.5*val*grad_gamma_diag1[ii,]
                // grad <- grad - 0.5*val*grad_gamma_diag1[jj,]
                for (int pp=0; pp<npar; pp++){
                    if ( grad[pp] != 0 ){
                        grad[pp] = grad[pp] / za;
                    }
                    if (grad_gamma_diag_bool(ii,pp)){
                        grad[pp] += -0.5*val*grad_gamma_diag1(ii,pp);
                    }
                    if (grad_gamma_diag_bool(jj,pp)){
                        grad[pp] += -0.5*val*grad_gamma_diag1(jj,pp);
                    }
                }

                // #-- discrepancy function
                x_ij = val;
                b0 = b0_jk(ii,jj);
                b1 = b1_jk(ii,jj);
                b2 = b2_jk(ii,jj);
                b3 = b3_jk(ii,jj);
                g1 = x_ij;
                g2 = x_ij*x_ij;
                g3 = g2*x_ij;

                // pm_exp <- b0.jk[ii,jj] + b1.jk[ii,jj]*x_ij + b2.jk[ii,jj]*x_ij^2 + b3.jk[ii,jj]*x_ij^3
                pm_exp = b0 + b1*g1 + b2*g2 + b3*g3;
                // pm_exp_der <- b1.jk[ii,jj] + 2*b2.jk[ii,jj]*x_ij + 3*b3.jk[ii,jj]*x_ij^2
                pm_exp_der = b1 + 2*b2*g1 + 3*b3*g2;
                // temp1 <- -2*wgtm[ii,jj] * ( pm[ii,jj] - pm_exp )
                temp1 = -2*wgtm(ii,jj)*(pm(ii,jj) - pm_exp)*pm_exp_der;
                for (int pp=0; pp<npar; pp++){
                    if (grad[pp] != 0){
                        grad[pp] = grad[pp]*temp1;
                    }
                }
                grad0 = grad0 + grad;
            } // end if
        } // end jj
    } // end ii


    //--- output
    return Rcpp::List::create(
            Rcpp::Named("grad_gamma_diag") = grad_gamma_diag,
            Rcpp::Named("grad0") = grad0
        );
}
///********************************************************************

