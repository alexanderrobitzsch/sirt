## File Name: noharm_sirt_optim_gradient_R_der_gamma_item.R
## File Version: 0.07


noharm_sirt_optim_gradient_R_der_gamma_item <- function(parm_table_free, Fmat,
    Pmat, FP, npar, NH, I)
{
    grad_gamma_diag <- matrix(0, nrow=I, ncol=npar)
    grad_gamma_diag_bool <- matrix(FALSE, nrow=I, ncol=npar)
    for (iii in 1L:I){
        grad <- rep(0,npar)
        grad_bool <- rep(FALSE,npar)
        for (hh in 1L:NH){
            mat_hh <- parm_table_free[hh,'mat']
            par_index_hh <- parm_table_free[hh,'index']
            row <- parm_table_free[hh,'row']
            col <- parm_table_free[hh,'col']
            # F
            if (mat_hh=='F'){
                if (row==iii){
                    der <- 2*FP[iii,col]
                    grad[par_index_hh] <- grad[par_index_hh] + der
                    grad_bool[par_index_hh] <- TRUE
                }
            }
            # P
            if (mat_hh=='P'){
                der <- 0
                if (row==col){
                    der <- Fmat[iii,col]^2
                    grad[par_index_hh] <- grad[par_index_hh] + der
                } else {
                    der <- Fmat[iii,row]*Fmat[iii,col]
                    grad[par_index_hh] <- grad[par_index_hh] + der
                }
                if (der!=0){ grad_bool[par_index_hh] <- TRUE }
            }
        }
        grad_gamma_diag[iii,] <- grad
        grad_gamma_diag_bool[iii,] <- grad_bool
    }

    #-- output
    res <- list(grad_gamma_diag=grad_gamma_diag,
                    grad_gamma_diag_bool=grad_gamma_diag_bool)
    return(res)
}
