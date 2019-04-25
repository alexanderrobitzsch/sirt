## File Name: lsem_weighted_mean.R
## File Version: 0.14

lsem_weighted_mean <- function( x, weights, x_resp=NULL )
{
    x <- as.matrix(x)
    if ( is.null(x_resp)){
        x_resp <- 1 - is.na(x)
    }
    weights_m <- weights * x_resp
    x[ ! x_resp ] <- 0
    weightsN <- colSums(weights_m)
    wm <- colSums( x * weights_m ) / ( weightsN - 1 )
    res <- list( weightsN=weightsN, mean=wm )
    return(res)
}
