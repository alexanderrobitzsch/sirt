## File Name: mgsem_scad_penalty.R
## File Version: 0.02
## File Last Change: 2022-02-25

mgsem_scad_penalty <- function(x, lambda, a=3.7)
{
    a <- max(a,1)
    x <- abs(x)
    res <- ifelse( x < lambda, lambda * x, 0)
    res <- res + ifelse( ( x >=lambda ) & ( x < a*lambda),
                        - ( x^2 - 2*a*lambda*x+lambda^2) / ( 2*(a-1)),0 )
    res <- res + ifelse( x>=a*lambda, (a+1)*lambda^2 / 2, 0 )
    return(res)
}
