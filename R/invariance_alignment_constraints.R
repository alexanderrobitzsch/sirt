## File Name: invariance_alignment_constraints.R
## File Version: 0.251
## File Last Change: 2019-03-06

invariance_alignment_constraints <- function(model, lambda_parm_tol, nu_parm_tol )
{
    CALL <- match.call()
    fpc_maxiter <- 20
    fpc_conv <- 1e-5

    nu <- model$nu.aligned
    lambda <- model$lambda.aligned
    wgt <- model$wgt
    miss_items <- model$miss_items

    #*** reparameterization
    lambda_list <- invariance_alignment_find_parameter_constraints(parm=lambda,
                        parm_tol=lambda_parm_tol, wgt=wgt, maxiter=fpc_maxiter,
                        conv=fpc_conv, miss_items=miss_items)
    nu_list <- invariance_alignment_find_parameter_constraints(parm=nu,
                        parm_tol=nu_parm_tol, wgt=wgt, maxiter=fpc_maxiter,
                        conv=fpc_conv, miss_items=miss_items)

    #--- output
    res <- list(model=model, nu_list=nu_list, lambda_list=lambda_list,
                    CALL=CALL)
    class(res) <- 'invariance_alignment_constraints'
    return(res)
}

