## File Name: dinvgamma2.R
## File Version: 0.01


dinvgamma2 <- function( x, n0, var0 )
{
    res <- sirt_dinvgamma( x, shape=n0 / 2, scale=n0 * var0 / 2 )
    return(res)
}
