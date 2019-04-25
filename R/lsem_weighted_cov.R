## File Name: lsem_weighted_cov.R
## File Version: 0.18

lsem_weighted_cov <- function( x, weights, x_resp=NULL )
{
    x <- as.matrix(x)
    if ( is.null(x_resp)){
        x_resp <- 1 - is.na(x)
    }
    eps0 <- 1e-100
    eps <- eps0 * max(weights)
    weights_m <- sqrt( weights + eps ) * x_resp
    x[ ! x_resp ] <- 0
    x_center <- lsem_weighted_mean( x=x, weights=weights_m, x_resp=x_resp)$mean
    XC <- matrix( x_center, nrow=nrow(x), ncol=ncol(x), byrow=TRUE )
    x <- x - XC
    weightsN <- crossprod(weights_m)
    xw <- as.matrix( x * weights_m)
    covw <- crossprod(xw) / weightsN
    Nobs <- mean( weightsN[ ! upper.tri(weightsN) ] )
    res <- list( weightsN=weightsN, cov=covw, Nobs=Nobs)
    return(res)
}
