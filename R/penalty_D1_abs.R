## File Name: penalty_D1_abs.R
## File Version: 0.02


penalty_D1_abs <- function(x, lambda, eps)
{
    res <- lambda*sqrt( x^2 + eps )
    return(res)
}
