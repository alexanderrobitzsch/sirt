## File Name: sirt_optimizer.R
## File Version: 0.304

sirt_optimizer <- function(optimizer, par, fn, grad=NULL, method="L-BFGS-B",
        hessian=TRUE, ...)
{
    s1 <- Sys.time()
    h <- 1e-5

    if (optimizer=="optim"){
        res <- stats::optim(par=par, fn=fn, gr=grad, method=method,
                        hessian=hessian, ...)
        res$iter <- res$counts['function']
    }
    if (optimizer=="nlminb"){
        res <- stats::nlminb(start=par, objective=fn, gradient=grad, ...)
        res$value <- res$objective
        a1 <- list(...)
        if (hessian){
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
        }
        res$iter <- res$iterations
    }
    s2 <- Sys.time()
    res$time_diff <- s2 - s1
    res$optimizer <- optimizer
    res$converged <- res$convergence==0
    #--- output
    return(res)
}
