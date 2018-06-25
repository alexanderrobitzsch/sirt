## File Name: sirt_optimizer.R
## File Version: 0.03

sirt_optimizer <- function(optimizer, par, fn, grad=NULL, method="L-BFGS-B", ...)
{
    s1 <- Sys.time()
    if (optimizer=="optim"){
        res <- stats::optim(par=par, fn=fn, grad=grad, method=method, ...)
    }
    if (optimizer=="nlminb"){
        res <- stats::nlminb(start=par, objective=fn, gradient=grad, ...)
    }
    s2 <- Sys.time()
    res$time_diff <- s2 - s1
    #--- output
    return(res)
}
