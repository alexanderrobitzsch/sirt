## File Name: mgsem_L0_penalty.R
## File Version: 0.03

mgsem_L0_penalty <- function(x=x, eps, gamma, deriv=FALSE, h=min(1e-4,eps/10) )
{
    if (deriv){
        y1 <- mgsem_L0_approx_ot(x=x+h, gamma=gamma, eps=eps)
        y2 <- mgsem_L0_approx_ot(x=x-h, gamma=gamma, eps=eps)
        y <- (y1-y2)/(2*h)
    } else {
        y <- mgsem_L0_approx_ot(x=x, gamma=gamma, eps=eps)
    }
    return(y)
}
