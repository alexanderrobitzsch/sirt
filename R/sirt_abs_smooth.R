## File Name: sirt_abs_smooth.R
## File Version: 0.07

sirt_abs_smooth <- function(x, deriv=0, eps=1e-4)
{
    .expr2 <- x^2 + eps
    y <- sqrt(.expr2)
    if (deriv==0){
        z <- y
    }
    if (deriv==1){
        z <- x / y
    }
    if (deriv==2){
        .expr4 <- 2 * x
        .expr5 <- 1 / y
        z <- .expr5 - x^2 * .expr2^(-1.5)
    }
    return(z)
}
