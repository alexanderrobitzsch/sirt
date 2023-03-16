## File Name: stratified_cronbach_alpha_compute_alpha.R
## File Version: 0.06
## File Last Change: 2020-02-10


#** compute alpha
stratified_cronbach_alpha_compute_alpha <- function( data )
{
    # covariance
    c1 <- stats::cov( data, use="pairwise.complete.obs" )
    # mean covariance
    c1a <- c1
    diag(c1a) <- 0
    I <- ncol(data)
    mc <- sum(c1a) / ( I^2 - I )

    # mean and variance
    mv <- mean( diag(c1) )
    alpha <- I*mc / ( mv + (I-1)*mc )
    mean.tot <- mean( rowSums(data), na.rm=TRUE )
    var.tot <- stats::var( rowSums(data), na.rm=TRUE )
    res <- list( I=I, alpha=alpha, mean.tot=mean.tot, var.tot=var.tot )
    return(res)
}
