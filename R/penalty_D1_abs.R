## File Name: penalty_D1_abs.R
## File Version: 0.02
## File Last Change: 2020-03-04


penalty_D1_abs <- function(x, lambda, eps)
{
    res <- lambda*sqrt( x^2 + eps )
    return(res)
}
