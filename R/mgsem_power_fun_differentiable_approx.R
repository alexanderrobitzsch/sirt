## File Name: mgsem_power_fun_differentiable_approx.R
## File Version: 0.096

mgsem_power_fun_differentiable_approx <- function(x, p, eps, deriv=FALSE,
        approx_method="lp")
{
    # logcomp <- TRUE
    logcomp <- FALSE
    if (deriv){  ## derivative
        if (approx_method=='lp'){
            if (p>0){
                if (!logcomp){
                    res <- p*((x^2+eps)^(p/2-1))*x
                } else {
                    logp <- log(p)
                    p2 <- (p/2-1)
                    res <- x*exp( p2*log(x^2+eps) + logp )
                }
            } else { # p=0
                res <- 2*x*eps/(x^2+eps)^2
            }
        }
        if (approx_method=='l2'){
            res <- 2*x
        }
    } else {  # no derivative
        if (approx_method=='lp'){
            if (p>0){
                res <- (x^2+eps)^(p/2)
            } else { # p=0
                res <- x^2 / (x^2 + eps)
            }
        }
        if (approx_method=='l2'){
            res <- x^2
        }
    }
    return(res)
}
