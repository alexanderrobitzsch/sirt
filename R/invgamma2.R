## File Name: invgamma2.R
## File Version: 0.207


#---- inverse gamma distribution for variance
rinvgamma2 <- function( n, n0, var0 )
{
    res <- sirt_rinvgamma( n, shape=n0 / 2, scale=n0 * var0 / 2 )
    return(res)
}



