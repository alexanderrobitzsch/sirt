## File Name: sirt_optimizer.R
## File Version: 0.13

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
        if (hessian){
            res$hessian <- CDM::numerical_gradient(par=res$par, FUN=grad, h=h)
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
