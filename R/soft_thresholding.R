## File Name: soft_thresholding.R
## File Version: 0.07


soft_thresholding <- function( x, lambda )
{
    x_abs <- abs(x)
    x <- ifelse( x_abs > lambda, x - sign(x) * lambda, 0 )
    return(x)
}


mlnormal_soft_thresholding <- soft_thresholding
