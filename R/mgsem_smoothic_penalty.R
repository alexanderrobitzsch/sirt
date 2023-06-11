## File Name: mgsem_smoothic_penalty.R
## File Version: 0.04


mgsem_smoothic_penalty <- function(x, eps, deriv=FALSE)
{
    if (deriv){  ## derivative
        res <- 2*x*eps / (x^2+eps)^2
        # 2*x*(x^2+eps) - x^2*(2*x)
    } else {  # no derivative
        res <- x^2 / (x^2 + eps)
    }
    return(res)
}
