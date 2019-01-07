## File Name: noharm_sirt_optim_function.R
## File Version: 0.18


noharm_sirt_optim_function <- function(x, parm_table, parm_index, I, D,
    b0.jk, b1.jk, b2.jk, b3.jk, pm, wgtm, use_rcpp=FALSE )
{

    parm_table <- noharm_sirt_partable_include_par(par=x, parm_table=parm_table)
    Fmat <- noharm_sirt_create_parameter_matrices("F", parm_table=parm_table, parm_index=parm_index)
    Pmat <- noharm_sirt_create_parameter_matrices("P", parm_table=parm_table, parm_index=parm_index)
    Psimat <- noharm_sirt_create_parameter_matrices("Psi", parm_table=parm_table, parm_index=parm_index)

    #- implied covariance: gamma values
    gamma_val <- noharm_sirt_implied_cov(Fmat=Fmat, Pmat=Pmat, Psimat=Psimat)

    ## gamma_ii
    gamma_diag <- diag(gamma_val)
    delta <- 1 + gamma_diag

    #- compute least squares function
    if (!use_rcpp){
        val <- noharm_sirt_optim_function_R( gamma_val=gamma_val, delta=delta, I=I,
                    wgtm=wgtm, pm=pm, b0.jk=b0.jk, b1.jk=b1.jk, b2.jk=b2.jk, b3.jk=b3.jk )
    } else {
        val <- sirt_rcpp_noharm_sirt_optim_fn_rcpp( gamma_val=gamma_val,
                    delta=delta, I=I, wgtm=wgtm, pm=pm, b0_jk=b0.jk, b1_jk=b1.jk, b2_jk=b2.jk,
                    b3_jk=b3.jk )
    }

    #-- output
    return(val)
}
