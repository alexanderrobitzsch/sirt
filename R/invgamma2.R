## File Name: invgamma2.R
## File Version: 0.13


###################################################################
# inverse gamma distribution for variance
rinvgamma2 <- function( n , n0 , var0 ){
    # INPUT:
    # N ... number of random draws
    # n0 ... sample size prior
    # var0 ... prior variance
#    res <- 1/ stats::rgamma( N , n0 / 2 ,  n0 * var0 / 2 )
    res <- sirt_rinvgamma( n , shape=n0 / 2 ,  scale=n0 * var0 / 2 )
    return(res)
}
#####################################################################
dinvgamma2 <- function( x , n0 , var0 ){
    res <- sirt_dinvgamma( x , shape=n0 / 2 ,  scale=n0 * var0 / 2 )
    return(res)
}

#-----------------------------------------------------------------------
# copied from MCMCpack package
sirt_dinvgamma <- function (x, shape, scale = 1)
{
    if (shape <= 0 | scale <= 0) {
        stop("Shape or scale parameter negative in dinvgamma().\n")
    }
    alpha <- shape
    beta <- scale
    log.density <- alpha * log(beta) - lgamma(alpha) - (alpha +
        1) * log(x) - (beta/x)
    return(exp(log.density))
}

sirt_rinvgamma <- function (n, shape, scale = 1)
{
    return(1/rgamma(n = n, shape = shape, rate = scale))
}
