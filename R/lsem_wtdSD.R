## File Name: lsem_wtdSD.R
## File Version: 0.16


lsem_wtdSD <- function( x, w )
{
    res1 <- sum( x*w )
    res2 <- sum( x^2*w)
    res12 <- res1^2
    if( res2 >=res12 ){
        res <- sqrt( res2 - res12 )
    } else {
        res <- 0
    }
    return(res)
}
