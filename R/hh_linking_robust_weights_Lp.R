## File Name: hh_linking_robust_weights_Lp.R
## File Version: 0.02


hh_linking_robust_weights_Lp <- function(x, p, eps)
{
    if (p==0){
        robw <- 2*eps/(x^2+eps)^2
    }
    if (p>0){
        robw <- p*(x^2+eps)^(p/2-1)
    }
    return(robw)
}
