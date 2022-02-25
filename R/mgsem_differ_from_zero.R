## File Name: mgsem_differ_from_zero.R
## File Version: 0.01

mgsem_differ_from_zero <- function(x, eps)
{
    res <- abs(x) > eps
    return(res)
}
