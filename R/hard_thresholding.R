

hard_thresholding <- function( x , lambda )
{
    x_abs <- abs(x)
    x <- ifelse( x_abs > lambda , x , 0 )
    return(x)
}
