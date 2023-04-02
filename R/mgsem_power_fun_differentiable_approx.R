## File Name: mgsem_power_fun_differentiable_approx.R
## File Version: 0.091

mgsem_power_fun_differentiable_approx <- function(x, p, eps, deriv=FALSE,
        approx_method="lp")
{
    # logcomp <- TRUE
    logcomp <- FALSE
    if (deriv){
        if (approx_method=='lp'){
            if (!logcomp){
                res <- p*((x^2+eps)^(p/2-1))*x
            } else {
                logp <- log(p)
                p2 <- (p/2-1)
                res <- x*exp( p2*log(x^2+eps) + logp )
            }
        }
        if (approx_method=='l2'){
            res <- 2*x
        }
    } else {
        if (approx_method=='lp'){
            res <- (x^2+eps)^(p/2)
        }
        if (approx_method=='l2'){
            res <- x^2
        }
    }
    return(res)
}
