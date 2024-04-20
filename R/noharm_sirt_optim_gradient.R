## File Name: noharm_sirt_optim_gradient.R
## File Version: 0.467


noharm_sirt_optim_gradient <- function(x, parm_table, parm_index, I, D,
        b0.jk, b1.jk, b2.jk, b3.jk, pm, wgtm, use_rcpp=FALSE )
{

    parm_table <- noharm_sirt_partable_include_par(par=x, parm_table=parm_table)
    Fmat <- noharm_sirt_create_parameter_matrices('F', parm_table=parm_table,
                            parm_index=parm_index)
    Pmat <- noharm_sirt_create_parameter_matrices('P', parm_table=parm_table,
                            parm_index=parm_index)
    Psimat <- noharm_sirt_create_parameter_matrices('Psi', parm_table=parm_table,
                            parm_index=parm_index)
    npar <- attr(parm_table, 'npar')

    # gamma values (diagonal)
    gamma_val <- noharm_sirt_implied_cov(Fmat=Fmat, Pmat=Pmat, Psimat=Psimat)
    FP <- Fmat %*% Pmat

    npar <- attr(parm_table, 'npar')
    NH <- attr(parm_table, 'NH')

    #- extract parameter table with free parameters
    parm_table_free <- parm_table[ attr(parm_table, 'parm_table_free_index'), ]

    #* computations
    if (use_rcpp){
        pt_index <- parm_table_free$index - 1
        pt_row <- parm_table_free$row -1
        pt_col <- parm_table_free$col - 1
        pt_matid <- parm_table_free$matid
        res <- sirt_rcpp_noharm_sirt_optim_gr_rcpp( gamma_val=gamma_val, NH=NH,
                    I=I, wgtm=wgtm, pm=pm, b0_jk=b0.jk, b1_jk=b1.jk, b2_jk=b2.jk,
                    b3_jk=b3.jk, npar=npar, pt_matid=pt_matid, pt_index=pt_index,
                    pt_row=pt_row, pt_col=pt_col, FP=FP, Fmat=Fmat, Pmat=Pmat,
                    Psimat=Psimat )
        grad0 <- res$grad0
    } else {
        grad0 <- noharm_sirt_optim_gradient_R( parm_table_free=parm_table_free,
                        Fmat=Fmat, Pmat=Pmat, Psimat=Psimat, FP=FP, npar=npar, NH=NH,
                        I=I, gamma_val=gamma_val, pm=pm, wgtm=wgtm, b0.jk=b0.jk,
                        b1.jk=b1.jk, b2.jk=b2.jk, b3.jk=b3.jk )
    }

    #-- output
    return(grad0)
}
