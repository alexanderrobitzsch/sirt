## File Name: rasch_mml2_difference_quotient.R
## File Version: 0.01

rasch_mml2_difference_quotient <- function(ll0, ll1, ll2, h)
{
    # first order derivative
    # f(x+h) - f(x-h) = 2* f'(x) * h    
    d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
    # second order derivative
    # f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
    d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
    #--- output
    res <- list( d1=d1, d2=d2)
    return(res)
}
