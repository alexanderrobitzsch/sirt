## File Name: noharm_sirt_optim_gradient_R_der_gamma_item_pair.R
## File Version: 0.18


noharm_sirt_optim_gradient_R_der_gamma_item_pair <- function(parm_table_free, Fmat,
    Pmat, Psimat, FP, npar, NH, I, gamma_val, grad_gamma_diag1, pm, b0.jk, b1.jk, b2.jk,
    b3.jk, wgtm, ii, jj )
{
    grad <- rep(0,npar)
    if (wgtm[ii,jj]>0){
        for (hh in 1:NH){
            mat_hh <- parm_table_free[hh,"mat"]
            par_index_hh <- parm_table_free[hh,"index"]
            row <- parm_table_free[hh,"row"]
            col <- parm_table_free[hh,"col"]
            # F
            if (mat_hh=="F"){
                der <- 0
                if (row==ii){
                    der <- FP[jj,col]
                    grad[par_index_hh] <- grad[par_index_hh] + der
                }
                if (row==jj){
                    der <- FP[ii,col]
                    grad[par_index_hh] <- grad[par_index_hh] + der
                }
            }
            # P
            if (mat_hh=="P"){
                der <- 0
                if (row==col){
                    der <- Fmat[ii,col]*Fmat[jj,col]
                    grad[par_index_hh] <- grad[par_index_hh] + der
                } else {
                    der <- Fmat[ii,row]*Fmat[jj,col]
                    grad[par_index_hh] <- grad[par_index_hh] + der
                }
            }
            # Psi
            if (mat_hh=="Psi"){
                if (row==ii){
                    if (col==jj){
                        der <- 1
                        grad[par_index_hh] <- grad[par_index_hh] + der
                    }
                }
            }
        }    # end hh
        # grad_gamma_nondiag <- grad

        # g_ij / sqrt( (1+g_ii)*(1+g_jj) )
        gii <- (1+gamma_val[ii,ii])
        gjj <- (1+gamma_val[jj,jj])
        gij <- gamma_val[ii,jj]
        za <- sqrt(gii*gjj)
        val <- gij / za

        # grad <- rep(0, npar)
        # first term
        grad <- grad / za
        # second term
        grad <- grad - 0.5*val*grad_gamma_diag1[ii,]
        # third term
        grad <- grad - 0.5*val*grad_gamma_diag1[jj,]

        #-- discrepancy function
        x_ij <- val
        pm_exp <- b0.jk[ii,jj] + b1.jk[ii,jj]*x_ij + b2.jk[ii,jj]*x_ij^2 + b3.jk[ii,jj]*x_ij^3
        pm_exp_der <- b1.jk[ii,jj] + 2*b2.jk[ii,jj]*x_ij + 3*b3.jk[ii,jj]*x_ij^2
        grad <- pm_exp_der*grad
        temp1 <- -2*wgtm[ii,jj] * ( pm[ii,jj] - pm_exp )
        grad <- temp1*grad
    }

    #--- output
    return(grad)
}
