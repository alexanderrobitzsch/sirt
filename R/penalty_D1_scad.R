## File Name: penalty_D1_scad.R
## File Version: 0.05


penalty_D1_scad <- function(x, lambda, eps, a=3.7)
{
    x <- abs(x)
    res <- ifelse( x < lambda, lambda * sqrt( x^2 + eps ), 0)
    res <- res + ifelse( ( x >=lambda ) & ( x < a*lambda),
                        - ( x^2 - 2*a*lambda*sqrt(x^2+eps)+lambda^2) / ( 2*(a-1)),0 )
    res <- res + ifelse( x>=a*lambda, (a+1)*lambda^2 / 2, 0 )
    return(res)
}
