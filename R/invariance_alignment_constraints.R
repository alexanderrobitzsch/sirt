## File Name: invariance_alignment_constraints.R
## File Version: 0.12

invariance_alignment_constraints <- function(model, lambda_parm_tol, nu_parm_tol )
{
    fpc_maxiter <- 20
    fpc_conv <- 1e-5

    nu <- model$nu.aligned
    lambda <- model$lambda.aligned
    wgt <- model$wgt

    #*** reparameterization
    lambda_list <- invariance_alignment_find_parameter_constraints(parm=lambda,
                        parm_tol=lambda_parm_tol, wgt=wgt, maxiter=fpc_maxiter, conv=fpc_conv)
    nu_list <- invariance_alignment_find_parameter_constraints(parm=nu,
                        parm_tol=nu_parm_tol, wgt=wgt, maxiter=fpc_maxiter, conv=fpc_conv)
    #--- output
    res <- list(model=model, nu_list=nu_list, lambda_list=lambda_list)
    class(res) <- 'invariance_alignment_constraints'
    return(res)
}

