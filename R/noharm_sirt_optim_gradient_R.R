## File Name: noharm_sirt_optim_gradient_R.R
## File Version: 0.03
## File Last Change: 2019-01-07

noharm_sirt_optim_gradient_R <- function(parm_table_free, Fmat, Pmat, Psimat, FP,
    npar, NH, I, gamma_val, pm, wgtm, b0.jk, b1.jk, b2.jk, b3.jk)
{

    #-- derivatives with respect to diagonale of gamma
    res <- noharm_sirt_optim_gradient_R_der_gamma_item( parm_table_free=parm_table_free,
                Fmat=Fmat, Pmat=Pmat, FP=FP, npar=npar, NH=NH, I=I )
    grad_gamma_diag <- res$grad_gamma_diag
    grad_gamma_diag_bool <- res$grad_gamma_diag_bool

    delta <- 1+diag(gamma_val)
    grad_gamma_diag1 <- grad_gamma_diag / delta
    grad0 <- rep(0, npar)
    for (ii in 1:(I-1)){
        for (jj in (ii+1):I){
            if (wgtm[ii,jj] >0 ){
                #- compute gradient for item pair
                grad1 <- noharm_sirt_optim_gradient_R_der_gamma_item_pair(
                            parm_table_free=parm_table_free, Fmat=Fmat, Pmat=Pmat,
                            Psimat=Psimat, FP=FP, npar=npar, NH=NH, I=I, gamma_val=gamma_val,
                            grad_gamma_diag1=grad_gamma_diag1, pm=pm, b0.jk=b0.jk,
                            b1.jk=b1.jk, b2.jk=b2.jk, b3.jk=b3.jk, wgtm=wgtm, ii=ii, jj=jj )
                grad0 <- grad0 + grad1
            }
        }
    }
    #--- output
    return(grad0)
}
