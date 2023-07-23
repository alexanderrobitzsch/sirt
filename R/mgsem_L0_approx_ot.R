## File Name: mgsem_L0_approx_ot.R
## File Version: 0.02


mgsem_L0_approx_ot <- function(x, gamma, eps)
{
    y <- 2/(1+exp(-gamma*(x^2+eps)^1/2) )-1
    return(y)
}
