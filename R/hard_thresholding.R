## File Name: hard_thresholding.R
## File Version: 0.06
## File Last Change: 2018-12-30


hard_thresholding <- function( x, lambda )
{
    x_abs <- abs(x)
    x <- ifelse( x_abs > lambda, x, 0 )
    return(x)
}
