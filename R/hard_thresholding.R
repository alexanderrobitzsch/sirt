## File Name: hard_thresholding.R
## File Version: 0.02
## File Last Change: 2017-01-18 11:02:47


hard_thresholding <- function( x , lambda )
{
    x_abs <- abs(x)
    x <- ifelse( x_abs > lambda , x , 0 )
    return(x)
}
