## File Name: sirt_optimizer_hessian.R
## File Version: 0.03

sirt_optimizer_hessian <- function(res, fn, grad, h)
{
    a1 <- list(...)
    arglist <- list(...)
    if (!is.null(grad)){
        fun <- grad
        hess_fun <- CDM::numerical_gradient
    } else {
        fun <- fn
        hess_fun <- CDM::numerical_Hessian
    }
    arglist <- sirt_remove_arguments_function(fun=fun, args=arglist)
    arglist$par <- res$par
    arglist$FUN <- fun
    arglist$h <- h
    res$hessian <- do.call(what=hess_fun, args=arglist)

    #--- output
    return(res)
}

