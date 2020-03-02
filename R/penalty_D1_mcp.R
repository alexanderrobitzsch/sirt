## File Name: penalty_D1_mcp.R
## File Version: 0.02


penalty_D1_mcp <- function(x, lambda, eps, a=2.7)
{
    x <- abs(x)
    res <- ifelse( x < a*lambda, lambda*sqrt(x^2 + eps) - x^2 / (2*a), .5*a*lambda^2)
    return(res)
}
